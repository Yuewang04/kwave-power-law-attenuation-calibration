function report = validate_kwave_cuda_y1_nearby(output_json)
%VALIDATE_KWAVE_CUDA_Y1_NEARBY Test exact and near-y=1 backend behavior.
%
% This isolated control extends validate_kwave_cuda_no_dispersion.m to
% y = 0.99, 1.00, and 1.01. Each MATLAB/CUDA and default/no_dispersion
% case is caught independently so an expected y=1 failure cannot stop the
% remaining runs. No production input or output is read or modified.

arguments
    output_json (1,1) string = ""
end

script_path = string(mfilename("fullpath"));
validation_dir = fileparts(script_path);
if strlength(output_json) == 0
    output_dir = fullfile(validation_dir, "results");
    if ~isfolder(output_dir), mkdir(output_dir); end
    stamp = string(datetime("now", "TimeZone", "UTC", "Format", "yyyyMMdd'T'HHmmss'Z'"));
    output_json = fullfile(output_dir, "kwave_cuda_y1_nearby_" + stamp + ".json");
else
    output_dir = fileparts(output_json);
    if strlength(output_dir) > 0 && ~isfolder(output_dir), mkdir(output_dir); end
end

kwave_root = string(getenv("TFUS_KWAVE_ROOT"));
if strlength(kwave_root) == 0
    error("Set TFUS_KWAVE_ROOT to the k-Wave Toolbox directory.");
end
addpath(genpath(kwave_root));
assert(exist("kspaceFirstOrder3D", "file") == 2, "MATLAB solver not found.");
assert(exist("kspaceFirstOrder3DG", "file") == 2, "CUDA wrapper not found.");

binary_path = string(getenv("TFUS_KWAVE_BINARY_PATH"));
if strlength(binary_path) == 0, binary_path = fullfile(kwave_root, "binaries"); end
binary_name = string(getenv("TFUS_KWAVE_BINARY_NAME"));
if strlength(binary_name) == 0
    if ispc, binary_name = "kspaceFirstOrder-CUDA.exe"; else, binary_name = "kspaceFirstOrder-CUDA"; end
end
binary_file = fullfile(binary_path, binary_name);

params = struct();
params.grid_size_xyz = [48, 20, 20];
params.spacing_m = 0.5e-3;
params.pml_xyz = [6, 6, 6];
params.cfl = 0.1;
params.duration_s = 18e-6;
params.frequency_hz = 500e3;
params.source_pressure_pa = 1e6;
params.sound_speed_m_s = 2800;
params.density_kg_m3 = 1900;
params.literature_alpha_coeff_db_cm_mhz = 1.86749;
params.exponents = [0.99, 1.00, 1.01];

report = struct();
report.schema_version = 1;
report.created_utc = char(datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
report.purpose = "Characterize exact and near-y=1 behavior in MATLAB and CUDA paths.";
report.scope = "Homogeneous small-grid control at y=0.99, 1.00, and 1.01; no production data used.";
report.matlab_version = version;
report.kwave_root = char(kwave_root);
report.kwave_version = read_kwave_version(kwave_root);
report.matlab_solver = which("kspaceFirstOrder3D");
report.cuda_wrapper = which("kspaceFirstOrder3DG");
report.cuda_binary = char(binary_file);
report.cuda_binary_exists = isfile(binary_file);
report.parameters = params;

[kgrid, source, sensor] = make_problem(params);
report.dt_s = kgrid.dt;
report.Nt = kgrid.Nt;

fields = struct();
case_reports = cell(numel(params.exponents), 1);
serialization_reports = cell(numel(params.exponents), 1);

for exponent_index = 1:numel(params.exponents)
    exponent = params.exponents(exponent_index);
    tag = exponent_tag(exponent);
    alpha_coeff = params.literature_alpha_coeff_db_cm_mhz * ...
        (params.frequency_hz / 1e6)^(1 - exponent);
    medium_default = make_medium(params, exponent, alpha_coeff, false);
    medium_no_disp = make_medium(params, exponent, alpha_coeff, true);

    fprintf("\n=== y=%.4f ===\n", exponent);
    [cpu_default_report, cpu_default_field] = run_matlab_case( ...
        "default", kgrid, medium_default, source, sensor, params);
    [cpu_no_disp_report, cpu_no_disp_field] = run_matlab_case( ...
        "no_dispersion", kgrid, medium_no_disp, source, sensor, params);

    if report.cuda_binary_exists
        [cuda_default_report, cuda_default_field] = run_cuda_case( ...
            "default", kgrid, medium_default, source, sensor, params, binary_path, binary_name);
        [cuda_no_disp_report, cuda_no_disp_field] = run_cuda_case( ...
            "no_dispersion", kgrid, medium_no_disp, source, sensor, params, binary_path, binary_name);
    else
        cuda_default_report = skipped_report("CUDA binary not found");
        cuda_no_disp_report = skipped_report("CUDA binary not found");
        cuda_default_field = [];
        cuda_no_disp_field = [];
    end

    fields.(tag).cpu_default = single(cpu_default_field);
    fields.(tag).cpu_no_dispersion = single(cpu_no_disp_field);
    fields.(tag).cuda_default = single(cuda_default_field);
    fields.(tag).cuda_no_dispersion = single(cuda_no_disp_field);

    one = struct();
    one.exponent = exponent;
    one.tan_pi_y_over_2 = tan(pi * exponent / 2);
    one.alpha_coeff_db_cm_mhzy = alpha_coeff;
    one.cpu_default = cpu_default_report;
    one.cpu_no_dispersion = cpu_no_disp_report;
    one.cuda_default = cuda_default_report;
    one.cuda_no_dispersion = cuda_no_disp_report;
    one.cpu_default_vs_no_dispersion = compare_if_available(cpu_default_field, cpu_no_disp_field);
    one.cuda_default_vs_no_dispersion = compare_if_available(cuda_default_field, cuda_no_disp_field);
    one.cpu_default_vs_cuda_default = compare_if_available(cpu_default_field, cuda_default_field);
    one.cpu_no_dispersion_vs_cuda_no_dispersion = compare_if_available(cpu_no_disp_field, cuda_no_disp_field);
    case_reports{exponent_index} = one;

    serialization_reports{exponent_index} = serialize_pair( ...
        tag, kgrid, medium_default, medium_no_disp, source, sensor, params, output_dir);
end

report.cases = case_reports;
report.serializations = serialization_reports;
report.interpretation_limits = [ ...
    "A nonzero finite result at y=0.99 or 1.01 does not establish physical " ...
    "accuracy. Exact-y=1 failure behavior is reported only for this tested " ...
    "MATLAB-wrapper/CUDA-binary pairing."];

mat_output = replace(output_json, ".json", "_pmax.mat");
save(mat_output, "fields", "params", "-v7.3");
report.pmax_mat_file = char(mat_output);
write_json(output_json, report);
fprintf("Wrote %s\n", output_json);
fprintf("Wrote %s\n", mat_output);
end

function [kgrid, source, sensor] = make_problem(params)
n = params.grid_size_xyz;
kgrid = kWaveGrid(n(1), params.spacing_m, n(2), params.spacing_m, n(3), params.spacing_m);
dt = params.cfl * params.spacing_m / params.sound_speed_m_s;
kgrid.setTime(ceil(params.duration_s / dt) + 1, dt);
source.p_mask = false(n);
source.p_mask(5, :, :) = true;
source.p = createCWSignals(kgrid.t_array, params.frequency_hz, params.source_pressure_pa, 0);
sensor.mask = true(n);
sensor.record = {'p_max'};
end

function medium = make_medium(params, exponent, alpha_coeff, no_dispersion)
medium.sound_speed = params.sound_speed_m_s;
medium.density = params.density_kg_m3;
medium.alpha_coeff = alpha_coeff;
medium.alpha_power = exponent;
if no_dispersion, medium.alpha_mode = 'no_dispersion'; end
end

function [result, field] = run_matlab_case(mode, kgrid, medium, source, sensor, params)
fprintf("Running MATLAB %s...\n", mode);
field = [];
started = tic;
try
    output = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', params.pml_xyz, ...
        'DataCast', 'single', 'PlotSim', false);
    runtime_s = toc(started);
    field = output.p_max;
    result = field_summary(field, runtime_s);
    result.status = "completed";
catch exception
    result = failure_report(exception, toc(started));
end
end

function [result, field] = run_cuda_case(mode, kgrid, medium, source, sensor, params, binary_path, binary_name)
fprintf("Running CUDA %s...\n", mode);
field = [];
started = tic;
try
    output = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', params.pml_xyz, ...
        'DataCast', 'single', 'BinaryPath', char(binary_path), ...
        'BinaryName', char(binary_name), 'DeleteData', true);
    runtime_s = toc(started);
    field = output.p_max;
    result = field_summary(field, runtime_s);
    result.status = "completed";
catch exception
    result = failure_report(exception, toc(started));
end
end

function result = field_summary(field, runtime_s)
values = double(field(:));
result = struct();
result.runtime_s = runtime_s;
result.peak_pressure_pa = max(values);
result.minimum_pressure_pa = min(values);
result.l2_norm = norm(values);
result.finite = all(isfinite(values));
result.nan_count = sum(isnan(values));
result.inf_count = sum(isinf(values));
end

function result = failure_report(exception, runtime_s)
result = struct();
result.status = "failed";
result.runtime_s = runtime_s;
result.error_identifier = exception.identifier;
result.error_message = exception.message;
end

function result = skipped_report(reason)
result = struct("status", "not_run", "reason", reason);
end

function result = compare_if_available(reference, candidate)
if isempty(reference) || isempty(candidate)
    result = struct("status", "not_comparable");
    return;
end
a = double(reference(:));
b = double(candidate(:));
delta = b - a;
scale = max(max(abs(a)), eps);
result = struct();
result.status = "completed";
result.exactly_equal = isequal(reference, candidate);
result.normalized_l2_difference = norm(delta) / max(norm(a), eps);
result.maximum_absolute_difference_pa = max(abs(delta));
result.maximum_relative_to_reference_peak = max(abs(delta)) / scale;
result.reference_peak_pressure_pa = max(a);
result.candidate_peak_pressure_pa = max(b);
result.peak_pressure_relative_difference = (max(b) - max(a)) / max(abs(max(a)), eps);
if std(a) > 0 && std(b) > 0
    correlation = corrcoef(a, b);
    result.spatial_correlation = correlation(1, 2);
else
    result.spatial_correlation = NaN;
end
end

function result = serialize_pair(tag, kgrid, medium_default, medium_no_disp, source, sensor, params, output_dir)
default_file = fullfile(output_dir, "serialization_" + tag + "_default.h5");
no_disp_file = fullfile(output_dir, "serialization_" + tag + "_no_dispersion.h5");
default_report = serialize_one(default_file, kgrid, medium_default, source, sensor, params);
no_disp_report = serialize_one(no_disp_file, kgrid, medium_no_disp, source, sensor, params);
result = struct();
result.exponent_tag = tag;
result.default = default_report;
result.no_dispersion = no_disp_report;
if default_report.status == "completed" && no_disp_report.status == "completed"
    result.comparison = compare_h5(default_file, no_disp_file);
else
    result.comparison = struct("status", "not_comparable");
end
end

function result = serialize_one(path, kgrid, medium, source, sensor, params)
if isfile(path), delete(path); end
try
    kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', params.pml_xyz, ...
        'DataCast', 'single', 'PlotSim', false, 'SaveToDisk', char(path));
    names = list_datasets(h5info(path), "");
    info = dir(path);
    result = struct();
    result.status = "completed";
    result.file = char(path);
    result.size_bytes = info.bytes;
    result.dataset_count = numel(names);
    result.dataset_names = cellstr(names);
catch exception
    result = failure_report(exception, 0);
    result.file = char(path);
    result.file_exists_after_failure = isfile(path);
end
end

function result = compare_h5(default_file, no_disp_file)
default_names = list_datasets(h5info(default_file), "");
no_disp_names = list_datasets(h5info(no_disp_file), "");
shared = intersect(default_names, no_disp_names, "stable");
different = strings(0, 1);
max_abs_diff = 0;
for index = 1:numel(shared)
    a = h5read(default_file, shared(index));
    b = h5read(no_disp_file, shared(index));
    if ~isequal(a, b)
        different(end + 1, 1) = shared(index); %#ok<AGROW>
        if isnumeric(a) && isnumeric(b) && isequal(size(a), size(b))
            max_abs_diff = max(max_abs_diff, max(abs(double(a(:)) - double(b(:))), [], "omitnan"));
        end
    end
end
result = struct();
result.status = "completed";
result.only_in_default = cellstr(setdiff(default_names, no_disp_names, "stable"));
result.only_in_no_dispersion = cellstr(setdiff(no_disp_names, default_names, "stable"));
result.different_datasets = cellstr(different);
result.maximum_absolute_dataset_difference = max_abs_diff;
result.dataset_payloads_identical = isempty(result.only_in_default) && ...
    isempty(result.only_in_no_dispersion) && isempty(result.different_datasets);
end

function names = list_datasets(group, prefix)
names = strings(0, 1);
for index = 1:numel(group.Datasets)
    names(end + 1, 1) = prefix + "/" + string(group.Datasets(index).Name); %#ok<AGROW>
end
for index = 1:numel(group.Groups)
    child = group.Groups(index);
    names = [names; list_datasets(child, string(child.Name))]; %#ok<AGROW>
end
names = unique(names, "stable");
end

function value = read_kwave_version(kwave_root)
value = "unknown";
contents_file = fullfile(kwave_root, "Contents.m");
if isfile(contents_file)
    source_text = fileread(contents_file);
    match = regexp(source_text, 'Version\s+([^\r\n]+)', 'tokens', 'once');
    if ~isempty(match), value = string(strtrim(match{1})); end
end
value = char(value);
end

function tag = exponent_tag(exponent)
tag = "y" + replace(compose("%.2f", exponent), ".", "p");
tag = char(tag);
end

function write_json(path, value)
json_text = jsonencode(value, "PrettyPrint", true);
fid = fopen(path, "w");
assert(fid >= 0, "Could not open JSON output: %s", path);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, json_text, "char");
end

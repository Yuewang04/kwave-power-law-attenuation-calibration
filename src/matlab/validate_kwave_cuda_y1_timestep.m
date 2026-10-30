function report = validate_kwave_cuda_y1_timestep(output_json)
%VALIDATE_KWAVE_CUDA_Y1_TIMESTEP Refine dt for default dispersion near y=1.
%
% Runs y=0.99 and y=1.01 with CFL 0.1, 0.05, and 0.025 on both the
% MATLAB and CUDA paths. This is an isolated homogeneous 3-D control and
% does not read or modify production data.

arguments
    output_json (1,1) string = ""
end

script_path = string(mfilename("fullpath"));
validation_dir = fileparts(script_path);
if strlength(output_json) == 0
    output_dir = fullfile(validation_dir, "results");
    if ~isfolder(output_dir), mkdir(output_dir); end
    stamp = string(datetime("now", "TimeZone", "UTC", "Format", "yyyyMMdd'T'HHmmss'Z'"));
    output_json = fullfile(output_dir, "kwave_cuda_y1_timestep_" + stamp + ".json");
else
    output_dir = fileparts(output_json);
    if strlength(output_dir) > 0 && ~isfolder(output_dir), mkdir(output_dir); end
end

kwave_root = string(getenv("TFUS_KWAVE_ROOT"));
if strlength(kwave_root) == 0, error("Set TFUS_KWAVE_ROOT."); end
addpath(genpath(kwave_root));
binary_path = string(getenv("TFUS_KWAVE_BINARY_PATH"));
binary_name = string(getenv("TFUS_KWAVE_BINARY_NAME"));
binary_file = fullfile(binary_path, binary_name);

params = struct();
params.grid_size_xyz = [48, 20, 20];
params.spacing_m = 0.5e-3;
params.pml_xyz = [6, 6, 6];
params.duration_s = 18e-6;
params.frequency_hz = 500e3;
params.source_pressure_pa = 1e6;
params.sound_speed_m_s = 2800;
params.density_kg_m3 = 1900;
params.literature_alpha_coeff_db_cm_mhz = 1.86749;
params.exponents = [0.99, 1.01];
params.cfl_values = [0.1, 0.05, 0.025];

report = struct();
report.schema_version = 1;
report.created_utc = char(datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
report.purpose = "Test timestep sensitivity of default dispersion at y=0.99 and y=1.01.";
report.scope = "Homogeneous small-grid control; MATLAB and CUDA; no production data used.";
report.matlab_version = version;
report.kwave_version = read_kwave_version(kwave_root);
report.cuda_binary = char(binary_file);
report.cuda_binary_exists = isfile(binary_file);
report.parameters = params;

fields = struct();
case_reports = cell(numel(params.exponents), numel(params.cfl_values));

for exponent_index = 1:numel(params.exponents)
    exponent = params.exponents(exponent_index);
    exponent_name = exponent_tag(exponent);
    alpha_coeff = params.literature_alpha_coeff_db_cm_mhz * ...
        (params.frequency_hz / 1e6)^(1 - exponent);
    medium.sound_speed = params.sound_speed_m_s;
    medium.density = params.density_kg_m3;
    medium.alpha_coeff = alpha_coeff;
    medium.alpha_power = exponent;

    for cfl_index = 1:numel(params.cfl_values)
        cfl = params.cfl_values(cfl_index);
        cfl_name = cfl_tag(cfl);
        [kgrid, source, sensor] = make_problem(params, cfl);
        fprintf("\n=== y=%.4f, CFL=%.4f, dt=%.9g s, Nt=%d ===\n", ...
            exponent, cfl, kgrid.dt, kgrid.Nt);

        [cpu_report, cpu_field] = run_matlab(kgrid, medium, source, sensor, params);
        if report.cuda_binary_exists
            [cuda_report, cuda_field] = run_cuda( ...
                kgrid, medium, source, sensor, params, binary_path, binary_name);
        else
            cuda_report = struct("status", "not_run", "reason", "CUDA binary not found");
            cuda_field = [];
        end

        fields.(exponent_name).(cfl_name).cpu_default = single(cpu_field);
        fields.(exponent_name).(cfl_name).cuda_default = single(cuda_field);
        one = struct();
        one.exponent = exponent;
        one.tan_pi_y_over_2 = tan(pi * exponent / 2);
        one.cfl = cfl;
        one.dt_s = kgrid.dt;
        one.Nt = kgrid.Nt;
        one.alpha_coeff_db_cm_mhzy = alpha_coeff;
        one.cpu_default = cpu_report;
        one.cuda_default = cuda_report;
        one.cpu_vs_cuda = compare_if_available(cpu_field, cuda_field);
        case_reports{exponent_index, cfl_index} = one;
    end
end

report.cases = case_reports;
report.refinement_summary = summarize_refinement(fields, params);
report.interpretation_limit = "Timestep behavior in this small control does not establish a universal safe exponent or CFL.";
mat_output = replace(output_json, ".json", "_pmax.mat");
save(mat_output, "fields", "params", "-v7.3");
report.pmax_mat_file = char(mat_output);
write_json(output_json, report);
fprintf("Wrote %s\n", output_json);
fprintf("Wrote %s\n", mat_output);
end

function [kgrid, source, sensor] = make_problem(params, cfl)
n = params.grid_size_xyz;
kgrid = kWaveGrid(n(1), params.spacing_m, n(2), params.spacing_m, n(3), params.spacing_m);
dt = cfl * params.spacing_m / params.sound_speed_m_s;
kgrid.setTime(ceil(params.duration_s / dt) + 1, dt);
source.p_mask = false(n);
source.p_mask(5, :, :) = true;
source.p = createCWSignals(kgrid.t_array, params.frequency_hz, params.source_pressure_pa, 0);
sensor.mask = true(n);
sensor.record = {'p_max'};
end

function [result, field] = run_matlab(kgrid, medium, source, sensor, params)
fprintf("Running MATLAB default...\n");
field = [];
started = tic;
try
    output = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', params.pml_xyz, ...
        'DataCast', 'single', 'PlotSim', false);
    field = output.p_max;
    result = field_summary(field, toc(started));
    result.status = "completed";
catch exception
    result = failure_report(exception, toc(started));
end
end

function [result, field] = run_cuda(kgrid, medium, source, sensor, params, binary_path, binary_name)
fprintf("Running CUDA default...\n");
field = [];
started = tic;
try
    output = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
        'PMLInside', false, 'PMLSize', params.pml_xyz, ...
        'DataCast', 'single', 'BinaryPath', char(binary_path), ...
        'BinaryName', char(binary_name), 'DeleteData', true);
    field = output.p_max;
    result = field_summary(field, toc(started));
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

function result = compare_if_available(reference, candidate)
if isempty(reference) || isempty(candidate)
    result = struct("status", "not_comparable");
    return;
end
a = double(reference(:));
b = double(candidate(:));
delta = b - a;
result = struct();
result.status = "completed";
result.exactly_equal = isequal(reference, candidate);
result.normalized_l2_difference = norm(delta) / max(norm(a), eps);
result.maximum_absolute_difference_pa = max(abs(delta));
result.peak_pressure_relative_difference = (max(b) - max(a)) / max(abs(max(a)), eps);
end

function summary = summarize_refinement(fields, params)
summary = cell(numel(params.exponents), 1);
for exponent_index = 1:numel(params.exponents)
    exponent = params.exponents(exponent_index);
    exponent_name = exponent_tag(exponent);
    cpu_peaks = NaN(size(params.cfl_values));
    cuda_peaks = NaN(size(params.cfl_values));
    for cfl_index = 1:numel(params.cfl_values)
        cfl_name = cfl_tag(params.cfl_values(cfl_index));
        cpu = fields.(exponent_name).(cfl_name).cpu_default;
        cuda = fields.(exponent_name).(cfl_name).cuda_default;
        if ~isempty(cpu), cpu_peaks(cfl_index) = max(double(cpu(:))); end
        if ~isempty(cuda), cuda_peaks(cfl_index) = max(double(cuda(:))); end
    end
    one = struct();
    one.exponent = exponent;
    one.cfl_values = params.cfl_values;
    one.cpu_peak_pressure_pa = cpu_peaks;
    one.cuda_peak_pressure_pa = cuda_peaks;
    one.cpu_peak_change_vs_previous = relative_changes(cpu_peaks);
    one.cuda_peak_change_vs_previous = relative_changes(cuda_peaks);
    summary{exponent_index} = one;
end
end

function changes = relative_changes(values)
changes = NaN(size(values));
for index = 2:numel(values)
    changes(index) = (values(index) - values(index - 1)) / max(abs(values(index - 1)), eps);
end
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

function tag = cfl_tag(cfl)
tag = "cfl" + replace(compose("%.3f", cfl), ".", "p");
tag = char(tag);
end

function write_json(path, value)
json_text = jsonencode(value, "PrettyPrint", true);
fid = fopen(path, "w");
assert(fid >= 0, "Could not open JSON output: %s", path);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, json_text, "char");
end

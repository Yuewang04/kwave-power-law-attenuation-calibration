function result = run_kwave_experiments(output_path, y, mode, coeff_scale, grid_factor, c_override, rho_override)
%RUN_KWAVE_EXPERIMENTS Homogeneous 1-D k-Wave PSTD attenuation benchmark.
%   Uses standard interior source placement (source_idx = 35*grid_factor,
%   outside the 20*grid_factor internal PML), full wavepacket duration
%   (180 us), exact Neper-to-dB conversion matching k-Wave's db2neper,
%   both positive and negative peak pressure tracking (p_max and p_min) to
%   catch non-oscillatory Hadamard blow-ups, phase-speed measurement, and
%   exact analytical PDE predictions (Treeby & Cox 2010, 2014).

if nargin < 1 || isempty(output_path)
    output_path = fullfile(pwd, 'results', 'matlab_results.json');
end
if nargin < 2 || isempty(y), y = 2; end
if nargin < 3, mode = ''; end
if nargin < 4 || isempty(coeff_scale), coeff_scale = 1; end
if nargin < 5 || isempty(grid_factor), grid_factor = 1; end
if nargin < 6 || isempty(c_override), c_override = 3100; end
if nargin < 7 || isempty(rho_override), rho_override = 2200; end

np2db = 20 * log10(exp(1));
f0 = 0.5;
alpha_lit = 21.5 * np2db / 100; % dB/(cm * MHz)
coeff_u = alpha_lit * coeff_scale;
coeff_c = alpha_lit * (f0^(1 - y)) * coeff_scale;

cal_u = run_one_kwave_case(f0, coeff_u, y, mode, grid_factor, c_override, rho_override);
cal_c = run_one_kwave_case(f0, coeff_c, y, mode, grid_factor, c_override, rho_override);

target_f0 = alpha_lit * f0 * coeff_scale;
cal_u = attach_target_metrics(cal_u, target_f0);
cal_c = attach_target_metrics(cal_c, target_f0);

freqs = [0.2 0.3 0.4 0.5 0.6 0.8 1.0];
sweep = [];
for i = 1:numel(freqs)
    item = run_one_kwave_case(freqs(i), coeff_c, y, mode, grid_factor, c_override, rho_override);
    item = attach_target_metrics(item, alpha_lit * freqs(i) * coeff_scale);
    if isempty(sweep)
        sweep = item;
    else
        sweep(end + 1) = item; %#ok<AGROW>
    end
end

result = struct();
result.constants = struct( ...
    'np2db_exact', np2db, ...
    'alpha_lit_Np_m_MHz', 21.5, ...
    'alpha_lit_dB_cm_MHz', alpha_lit, ...
    'f0_MHz', f0, ...
    'target_dB_cm_at_f0', target_f0, ...
    'alpha_coeff_uncorrected', coeff_u, ...
    'alpha_coeff_corrected', coeff_c, ...
    'alpha_power', y, ...
    'coeff_scale', coeff_scale, ...
    'grid_factor', grid_factor, ...
    'sound_speed_m_s', c_override, ...
    'density_kg_m3', rho_override, ...
    'source_idx_base', 35, ...
    'pml_size_base', 20, ...
    't_end_s', 180e-6);
result.calibration = struct('uncorrected', cal_u, 'corrected', cal_c);
result.frequency_sweep = sweep;

out_dir = fileparts(output_path);
if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
fid = fopen(output_path, 'w');
assert(fid >= 0, 'Cannot open output path: %s', output_path);
cleanup_obj = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(result, 'PrettyPrint', true), 'char');
end

function item = attach_target_metrics(item, target_dB_cm)
item.target_dB_cm = target_dB_cm;
if strcmp(item.field_status, 'stable') && isfinite(item.attenuation_dB_cm)
    item.relative_error_percent = 100 * (item.attenuation_dB_cm - target_dB_cm) / target_dB_cm;
else
    item.relative_error_percent = NaN;
end
if isfinite(item.analytical_pde_attenuation_dB_cm)
    item.analytical_pde_error_percent = 100 * (item.analytical_pde_attenuation_dB_cm - target_dB_cm) / target_dB_cm;
else
    item.analytical_pde_error_percent = NaN;
end
end

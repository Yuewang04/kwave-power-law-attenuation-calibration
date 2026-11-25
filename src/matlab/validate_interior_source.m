function result = validate_interior_source(output_path)
% VALIDATE_INTERIOR_SOURCE
% Standard k-Wave 1-D benchmark with the source placed strictly in the
% physical interior (grid index 35, well outside the 20-point internal PML),
% and recording time extended to 150 us to ensure full waveform passage.
%
% This script validates that the ~50% parameter mismatch at y=2 and the ~57%
% dispersion deficit near y=1 are intrinsic properties of k-Wave rather than
% artifacts of placing the source near/inside the PML.

if nargin < 1
    output_path = fullfile(pwd, 'results', 'interior_source_validation.json');
end

f0 = 0.5;
np2db = 20 * log10(exp(1));
alpha_lit = 21.5 * np2db / 100; % dB/(cm MHz)
target_f0 = alpha_lit * f0;

fprintf('=== Standard Interior-Source k-Wave Validation ===\n');
fprintf('Source index: 35 (PML size: 20 grid points)\n');
fprintf('Recording time: 150 us (full wave packet transmission)\n');
fprintf('Literature target at 0.5 MHz: %.4f dB/cm\n\n', target_f0);

% 1. y=2 uncorrected
coeff_u = alpha_lit;
r_u = run_single_interior(f0, coeff_u, 2, '');
err_u = 100 * (r_u.attenuation_dB_cm - target_f0) / target_f0;
fprintf('1. y=2 Uncorrected: fitted = %.4f dB/cm (target = %.4f, error = %+.2f%%)\n', ...
    r_u.attenuation_dB_cm, target_f0, err_u);

% 2. y=2 corrected
coeff_c = alpha_lit * f0^(1 - 2);
r_c = run_single_interior(f0, coeff_c, 2, '');
err_c = 100 * (r_c.attenuation_dB_cm - target_f0) / target_f0;
fprintf('2. y=2 Corrected:   fitted = %.4f dB/cm (target = %.4f, error = %+.2f%%)\n', ...
    r_c.attenuation_dB_cm, target_f0, err_c);

% 3. y=1.01 near singularity
coeff_101 = alpha_lit * f0^(1 - 1.01);
r_101 = run_single_interior(f0, coeff_101, 1.01, '');
err_101 = 100 * (r_101.attenuation_dB_cm - target_f0) / target_f0;
fprintf('3. y=1.01 Singular: fitted = %.4f dB/cm (target = %.4f, error = %+.2f%%)\n\n', ...
    r_101.attenuation_dB_cm, target_f0, err_101);

result = struct(...
    'description', 'Standard interior-source control (source_idx=35, pml_size=20, t_end=150us)', ...
    'f0_MHz', f0, ...
    'target_dB_cm', target_f0, ...
    'uncorrected_y2', struct('fitted_dB_cm', r_u.attenuation_dB_cm, 'error_percent', err_u), ...
    'corrected_y2', struct('fitted_dB_cm', r_c.attenuation_dB_cm, 'error_percent', err_c), ...
    'near_y1p01', struct('fitted_dB_cm', r_101.attenuation_dB_cm, 'error_percent', err_101) ...
);

fid = fopen(output_path, 'w');
if fid ~= -1
    fwrite(fid, jsonencode(result), 'char');
    fclose(fid);
    fprintf('Results saved to %s\n', output_path);
end
end

function r = run_single_interior(f_MHz, coeff, y, mode_override)
Nx = 512; dx = 0.5e-3; c = 3100; rho = 2200; f = f_MHz * 1e6;
kgrid = kWaveGrid(Nx, dx);
kgrid.makeTime(c, [], 180e-6);

medium.sound_speed = c;
medium.density = rho;
medium.alpha_coeff = coeff;
medium.alpha_power = y;
if y == 1
    medium.alpha_mode = 'no_dispersion';
elseif ~isempty(mode_override)
    medium.alpha_mode = mode_override;
end

source_idx = 35;
source.p_mask = zeros(Nx, 1);
source.p_mask(source_idx) = 1;
source.p = toneBurst(1/kgrid.dt, f, 18, 'SignalOffset', 0);

pml_size = 20;
rel_offsets = [52, 92, 132, 172, 212, 252]';
idx = source_idx + rel_offsets;
sensor.mask = false(Nx, 1);
sensor.mask(idx) = true;
sensor.record = {'p_max'};

data = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...
    'PMLSize', pml_size, 'PlotSim', false, 'DataCast', 'double', 'Smooth', false);

amp = abs(double(data.p_max(:)));
x = (idx - source_idx) * dx;
good = amp > 0;
fit = polyfit(x(good), log(amp(good)), 1);
slope = fit(1);
attenuation_dB_cm = (-slope) * (20 * log10(exp(1))) / 100;
r = struct('attenuation_dB_cm', attenuation_dB_cm);
end

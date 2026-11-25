function r = run_one_kwave_case(f_MHz, coeff, y, mode_override, grid_factor, c, rho, dt_factor)
%RUN_ONE_KWAVE_CASE Single 1-D k-Wave simulation with full stability & phase diagnostics.
%   Places source at index 35*grid_factor (15*grid_factor cells inside physical
%   interior beyond the 20*grid_factor internal PML), runs for 180 us so even
%   dispersively slowed 0.2-MHz 18-cycle tonebursts clear the farthest sensor
%   (126 mm), records {'p', 'p_max', 'p_min'}, detects Hadamard blow-up via
%   max(|p|), measures unaliased phase speed using an adjacent +dx pair, and
%   computes exact analytical fractional-Laplacian PDE predictions.

if nargin < 4, mode_override = ''; end
if nargin < 5 || isempty(grid_factor), grid_factor = 1; end
if nargin < 6 || isempty(c), c = 3100; end
if nargin < 7 || isempty(rho), rho = 2200; end
if nargin < 8 || isempty(dt_factor), dt_factor = 1; end

np2db = 20 * log10(exp(1));
Nx = 512 * grid_factor;
dx = 0.5e-3 / grid_factor;
f = f_MHz * 1e6;
t_end = 180e-6;

kgrid = kWaveGrid(Nx, dx);
kgrid.makeTime(c, [], t_end);
if dt_factor ~= 1
    dt_base = kgrid.dt;
    Nt_base = kgrid.Nt;
    kgrid.setTime(round(Nt_base * dt_factor), dt_base / dt_factor);
end

medium = struct();
medium.sound_speed = c;
medium.density = rho;
medium.alpha_coeff = coeff;
medium.alpha_power = y;

effective_mode = 'default_dispersive';
if y == 1
    medium.alpha_mode = 'no_dispersion';
    effective_mode = 'no_dispersion';
elseif ~isempty(mode_override)
    medium.alpha_mode = mode_override;
    effective_mode = mode_override;
end

pml_size = 20 * grid_factor;
source_idx = 35 * grid_factor;
rel_offsets = [52, 92, 132, 172, 212, 252]' * grid_factor;
idx = source_idx + rel_offsets;
% Add adjacent sensor at idx(1) + grid_factor (spacing = 0.5 mm << lambda/2)
% to measure unaliased phase speed cp(f) cleanly.
phase_pair_idx = idx(1) + grid_factor;
all_idx = [idx(1); phase_pair_idx; idx(2:end)];

source = struct();
source.p_mask = zeros(Nx, 1);
source.p_mask(source_idx) = 1;
source.p = toneBurst(1 / kgrid.dt, f, 18, 'SignalOffset', 0);

sensor = struct();
sensor.mask = false(Nx, 1);
sensor.mask(all_idx) = true;
sensor.record = {'p', 'p_max', 'p_min'};

data = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...
    'PMLSize', pml_size, 'PlotSim', false, 'DataCast', 'double', 'Smooth', false);

p_pos_all = double(data.p_max(:));
p_neg_all = double(data.p_min(:));
% Extract the 6 attenuation sensors: indices [1, 3, 4, 5, 6, 7]
att_sel = [1, 3, 4, 5, 6, 7];
p_pos = p_pos_all(att_sel);
p_neg = p_neg_all(att_sel);
amp_abs = max(abs(p_pos), abs(p_neg));
max_abs_p = max(max(abs(p_pos_all)), max(abs(p_neg_all)));

wave_all = double(data.p);
if size(wave_all, 1) == numel(all_idx)
    wave_all = wave_all.'; % Nt x 7
end
wave = wave_all(:, att_sel);

x = (idx - source_idx) * dx;

% Exact analytical PDE operator quantities (Treeby & Cox 2010, 2014)
a0_np = db2neper(coeff, y);
w = 2 * pi * f;
k0 = w / c;
if strcmp(effective_mode, 'no_dispersion')
    eta = 0;
    stiffness_M_k0 = 1.0;
    treeby_denom = 1.0;
    analytical_att_dB_cm = coeff * (f_MHz ^ y);
    analytical_cp_m_s = c;
else
    tan_val = tan(pi * y / 2);
    eta = 2 * a0_np * (c ^ y) * tan_val;
    stiffness_M_k0 = 1 - eta * (k0 ^ (y - 1));
    treeby_denom = 1 - (y + 1) * a0_np * c * tan_val * (w ^ (y - 1));
    if stiffness_M_k0 > 0 && treeby_denom > 0
        analytical_att_dB_cm = (a0_np * (w ^ y) / treeby_denom) * (0.01 * np2db);
        analytical_cp_m_s = c * sqrt(stiffness_M_k0);
    else
        analytical_att_dB_cm = NaN;
        analytical_cp_m_s = NaN;
    end
end

% Check stability: a passive 1-Pa toneburst in an absorbing medium must stay O(1) Pa.
% Exponential Hadamard blow-up (when stiffness_M_k0 <= 0) drives |p| >> 1e6 Pa or Inf/NaN.
is_diverged = any(~isfinite(amp_abs)) || (max_abs_p > 1e6) || any(~isfinite(wave_all(:)));

if is_diverged
    field_status = 'diverged';
    attenuation_dB_cm = NaN;
    slope = NaN;
    intercept = NaN;
    r2 = NaN;
    ln_p = NaN(1, numel(idx));
    fit_res = NaN(1, numel(idx));
    att_spec_dB_cm = NaN;
    r2_spec = NaN;
    cp_measured = NaN;
    spec_freq_MHz = NaN;
else
    field_status = 'stable';
    good = amp_abs > 0;
    ln_col = log(amp_abs(good));
    fit_v = polyfit(x(good), ln_col, 1);
    slope = fit_v(1);
    intercept = fit_v(2);
    attenuation_dB_cm = (-slope) * np2db / 100;
    pred = slope * x(good) + intercept;
    ss = sum((ln_col - mean(ln_col)).^2);
    sr = sum((ln_col - pred).^2);
    r2 = 1 - sr / max(ss, realmin);
    ln_p = log(amp_abs)';
    fit_res = (ln_col - pred)';

    % Spectral amplitude across 6 sensors & unaliased phase speed across pair (1, 2)
    Nt_actual = size(wave_all, 1);
    freq_axis = (0:Nt_actual-1)' / (Nt_actual * kgrid.dt);
    half_n = floor(Nt_actual / 2) + 1;
    [~, bin] = min(abs(freq_axis(1:half_n) - f));
    spec_freq_MHz = freq_axis(bin) / 1e6;
    fft_all = fft(wave_all, [], 1);
    spec_amp = abs(fft_all(bin, att_sel))';

    if all(spec_amp > 0)
        ln_spec = log(spec_amp);
        fit_s = polyfit(x, ln_spec, 1);
        att_spec_dB_cm = (-fit_s(1)) * np2db / 100;
        pred_s = fit_s(1) * x + fit_s(2);
        ss_s = sum((ln_spec - mean(ln_spec)).^2);
        sr_s = sum((ln_spec - pred_s).^2);
        r2_spec = 1 - sr_s / max(ss_s, realmin);
    else
        att_spec_dB_cm = NaN;
        r2_spec = NaN;
    end

    % Unaliased phase difference across 0.5 mm (spacing = grid_factor * dx)
    dx_pair = (phase_pair_idx - idx(1)) * dx;
    dphi = abs(angle(fft_all(bin, 2) / fft_all(bin, 1)));
    cp_measured = (2 * pi * freq_axis(bin) * dx_pair) / max(dphi, realmin);
end

r = struct( ...
    'frequency_MHz', f_MHz, ...
    'alpha_coeff', coeff, ...
    'alpha_power', y, ...
    'alpha_mode', effective_mode, ...
    'cycles', 18, ...
    'field_status', field_status, ...
    'max_abs_pressure_Pa', sanitize_nonfinite(max_abs_p), ...
    'attenuation_dB_cm', attenuation_dB_cm, ...
    'attenuation_spectral_dB_cm', att_spec_dB_cm, ...
    'fit_R2', r2, ...
    'fit_R2_spectral', r2_spec, ...
    'spectral_frequency_MHz', spec_freq_MHz, ...
    'measured_phase_speed_m_s', cp_measured, ...
    'analytical_pde_attenuation_dB_cm', analytical_att_dB_cm, ...
    'analytical_phase_speed_m_s', analytical_cp_m_s, ...
    'absorb_eta', eta, ...
    'stiffness_M_k0', stiffness_M_k0, ...
    'treeby_denominator', treeby_denom, ...
    'distances_m', x', ...
    'p_max_pos', sanitize_vec(p_pos'), ...
    'p_min_neg', sanitize_vec(p_neg'), ...
    'p_max', sanitize_vec(amp_abs'), ...
    'ln_p_max', ln_p, ...
    'fit_residual', fit_res, ...
    'fit_slope_Np_per_m', slope, ...
    'fit_intercept', intercept, ...
    'source_idx_grid', source_idx, ...
    'pml_size_grid', pml_size, ...
    'domain_length_m', Nx * dx, ...
    'dt_s', kgrid.dt, ...
    'Nt', kgrid.Nt, ...
    't_end_s', t_end, ...
    'grid_factor', grid_factor, ...
    'wavelength_m', c / f);
end

function v = sanitize_nonfinite(x)
if isinf(x)
    v = sign(x) * 1.7976931348623157e308;
else
    v = x;
end
end

function v = sanitize_vec(arr)
v = arr;
inf_mask = isinf(v);
if any(inf_mask)
    v(inf_mask) = sign(v(inf_mask)) * 1.7976931348623157e308;
end
end

function run_kwave_experiments(output_path, y, mode, coeff_scale, grid_factor)
% Real MATLAB k-Wave 1-D PSTD experiments. alpha_coeff is dB/(MHz^y cm).
if nargin < 1, output_path = fullfile(pwd, 'results', 'matlab_results.json'); end
if nargin < 2, y = 2; end
% Optional third argument lets diagnostics force the k-Wave absorption mode.
% Leave empty for the normal power-law dispersive implementation.
mode_override = '';
if nargin >= 3, mode_override = mode; end
if nargin < 4, coeff_scale = 1; end
if nargin < 5, grid_factor = 1; end
f0 = 0.5; alpha_lit = 21.5 * 8.686 / 100; coeff_u = alpha_lit; coeff_c = alpha_lit * f0^(1-y);
coeff_u = coeff_u * coeff_scale; coeff_c = coeff_c * coeff_scale;
cal_u = run_one(f0, coeff_u, y, mode_override, grid_factor); cal_c = run_one(f0, coeff_c, y, mode_override, grid_factor);
freqs = [0.2 0.3 0.4 0.5 0.6 0.8 1.0]; sweep = [];
for i = 1:numel(freqs)
    item = run_one(freqs(i), coeff_c, y, mode_override, grid_factor);
    item.target_dB_cm = alpha_lit * freqs(i);
    item.relative_error_percent = 100 * (item.attenuation_dB_cm - item.target_dB_cm) / item.target_dB_cm;
    if isempty(sweep), sweep = item; else, sweep(end+1) = item; end
end
cal_u.target_dB_cm = alpha_lit*f0; cal_c.target_dB_cm = alpha_lit*f0;
cal_u.relative_error_percent = 100*(cal_u.attenuation_dB_cm-cal_u.target_dB_cm)/cal_u.target_dB_cm;
cal_c.relative_error_percent = 100*(cal_c.attenuation_dB_cm-cal_c.target_dB_cm)/cal_c.target_dB_cm;
result.constants = struct('alpha_lit_dB_cm_MHz',alpha_lit,'f0_MHz',f0,'target_dB_cm_at_f0',alpha_lit*f0,'alpha_coeff_uncorrected',coeff_u,'alpha_coeff_corrected',coeff_c,'alpha_power',y,'coeff_scale',coeff_scale,'grid_factor',grid_factor);
result.calibration = struct('uncorrected',cal_u,'corrected',cal_c); result.frequency_sweep = sweep;
fid=fopen(output_path,'w'); fwrite(fid,jsonencode(result),'char'); fclose(fid);
end

function r = run_one(f_MHz, coeff, y, mode_override, grid_factor)
Nx=512*grid_factor; dx=0.5e-3/grid_factor; c=3100; rho=2200; f=f_MHz*1e6;
kgrid=kWaveGrid(Nx,dx); kgrid.makeTime(c, [], 90e-6);
medium.sound_speed=c; medium.density=rho; medium.alpha_coeff=coeff; medium.alpha_power=y;
if y == 1
    medium.alpha_mode = 'no_dispersion';
elseif ~isempty(mode_override)
    medium.alpha_mode = mode_override;
end
source_idx=8*grid_factor; pml_size=20*grid_factor;
% NOTE: with pml_size=20 the source at index 8 sits *inside* the internal PML.
% This is not standard k-Wave practice and is disclosed in the README Limitations.
% Moving the source into the physical interior is a required follow-up before these
% numbers are used for anything quantitative.
source.p_mask=zeros(Nx,1); source.p_mask(source_idx)=1; source.p=toneBurst(1/kgrid.dt,f,18,'SignalOffset',0);
sensor.mask=false(Nx,1); idx=(60:40:280)'*grid_factor; sensor.mask(idx)=true; sensor.record={'p_max'};
data=kspaceFirstOrder1D(kgrid,medium,source,sensor,'PMLSize',pml_size,'PlotSim',false,'DataCast','double','Smooth',false);
amp=abs(double(data.p_max(:))); x=(idx-source_idx)*dx; good=amp>0;
fit=polyfit(x(good),log(amp(good)),1); slope=fit(1); attenuation_dB_cm=(-slope)*8.686/100;
pred=fit(1)*x(good)+fit(2); ss=sum((log(amp(good))-mean(log(amp(good)))).^2); sr=sum((log(amp(good))-pred).^2); r2=1-sr/ss;
r=struct('frequency_MHz',f_MHz,'alpha_coeff',coeff,'alpha_power',y,'cycles',18,'attenuation_dB_cm',attenuation_dB_cm,'distances_m',x','p_max',amp','ln_p_max',log(amp)','fit_residual',log(amp(good))'-pred','fit_slope_Np_per_m',slope,'fit_intercept',fit(2),'fit_R2',r2,'domain_length_m',Nx*dx,'pml_size_grid',pml_size,'grid_factor',grid_factor,'wavelength_m',c/f);
end

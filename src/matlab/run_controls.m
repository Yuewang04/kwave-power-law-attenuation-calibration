function run_controls(output_path)
% Focused controls for the y=0.98 peak-pressure sensitivity observation.
if nargin < 1, output_path = fullfile(pwd, 'results', 'controls_y098.json'); end
f0 = 0.5; y = 0.98; c = 3100; rho = 2200;
alpha_lit = 21.5 * 8.686 / 100; coeff = alpha_lit * f0^(1-y);
dt_factors = [1 2 4]; controls_cell = cell(1,numel(dt_factors));
for k = 1:numel(dt_factors)
    controls_cell{k} = run_case(f0, coeff, y, c, rho, alpha_lit, dt_factors(k));
end
controls = [controls_cell{:}];
result = struct('case','homogeneous_1D_y0p98_corrected_default_dispersion', ...
    'constants',struct('f0_MHz',f0,'alpha_lit_dB_cm_MHz',alpha_lit, ...
    'target_dB_cm_at_f0',alpha_lit*f0,'alpha_coeff',coeff,'alpha_power',y), ...
    'controls',controls);
fid=fopen(output_path,'w'); fwrite(fid,jsonencode(result),'char'); fclose(fid);
end

function r = run_case(f_MHz, coeff, y, c, rho, alpha_lit, dt_factor)
Nx=512; dx=0.5e-3; f=f_MHz*1e6;
kgrid=kWaveGrid(Nx,dx); kgrid.makeTime(c, [], 90e-6);
dt_base=kgrid.dt; Nt_base=kgrid.Nt;
kgrid.setTime(round(Nt_base*dt_factor), dt_base/dt_factor);
medium.sound_speed=c; medium.density=rho; medium.alpha_coeff=coeff; medium.alpha_power=y;
source_idx=8; p_mask=zeros(Nx,1); p_mask(source_idx)=1;
source.p_mask=p_mask; source.p=toneBurst(1/kgrid.dt,f,18,'SignalOffset',0);
sensor.mask=false(Nx,1); idx=(60:40:280)'; sensor.mask(idx)=true; sensor.record={'p'};
data=kspaceFirstOrder1D(kgrid,medium,source,sensor,'PMLSize',20,'PlotSim',false,'DataCast','double','Smooth',false);
wave=double(data.p);
if size(wave,1) == numel(idx), wave=wave.'; end
x=((idx-source_idx)*dx); x=x(:);
% Match the existing estimator: positive peak returned by p_max, then abs().
amp=max(wave,[],1); amp=abs(amp(:));
[att_peak,r2_peak,slope,intc]=fit_attenuation(x,amp);
% Frequency-domain control: amplitude at the nearest FFT bin to f0.
freq_axis=(0:size(wave,1)-1)'/(size(wave,1)*kgrid.dt);
[~,bin]=min(abs(freq_axis-f)); spec=zeros(numel(idx),1);
for j=1:numel(idx)
    fft_wave=fft(wave(:,j),[],1);
    spec(j)=abs(fft_wave(bin));
end
[att_spec,r2_spec,~,~]=fit_attenuation(x,spec);
r=struct('dt_factor',dt_factor,'dt_s',kgrid.dt,'Nt',kgrid.Nt, ...
 'attenuation_peak_dB_cm',att_peak,'relative_error_peak_percent',100*(att_peak-alpha_lit*f_MHz)/ (alpha_lit*f_MHz), ...
 'R2_peak',r2_peak,'attenuation_spectral_dB_cm',att_spec,'relative_error_spectral_percent',100*(att_spec-alpha_lit*f_MHz)/(alpha_lit*f_MHz), ...
 'R2_spectral',r2_spec,'spectral_frequency_MHz',freq_axis(bin)/1e6,'distances_m',x');
end

function [att,r2,slope,intercept]=fit_attenuation(x,amp)
good=amp>0; fitv=polyfit(x(good),log(amp(good)),1); slope=fitv(1); intercept=fitv(2);
att=(-slope)*8.686/100; pred=fitv(1)*x(good)+fitv(2);
ss=sum((log(amp(good))-mean(log(amp(good)))).^2); sr=sum((log(amp(good))-pred).^2); r2=1-sr/ss;
end

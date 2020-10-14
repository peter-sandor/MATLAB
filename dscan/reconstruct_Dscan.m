%% Load Fundamental Spectrum (own measurement)
filename_sp='U:\Measurement_Data\OAC_IFROG\2020_10_06_ACER\Venteon_spectrum.txt';

if 1
    data_sp = read_spectrum(filename_sp);
    data_sp_norm = condition_spectrum(data_sp);
    wl = data_sp_norm(:,1);
    I_wl = data_sp_norm(:,2);
else
    venteon_spectrum=load(filename_sp);
    wl0=venteon_spectrum(:,1);
    sp0=venteon_spectrum(:,2);
    sp0=(sp0-mean(sp0(1:50)))/max(sp0-mean(sp0(1:50)));
    crop_start=1; crop_end=1024;
    sp0(sp0<0.01)=0;
    wl=map2colvec(wl0(crop_start:crop_end));
    I_wl=map2colvec(sp0(crop_start:crop_end));
end
phi_wl = zeros([length(wl) 1]);

E_om = SpectrumConvert([map2colvec(wl) map2colvec(I_wl) map2colvec(phi_wl)]);
[omega_center,sigma_Venteon] = peak_props([E_om(:,1) abs(E_om(:,2)).^2]);
E_om(:,3) = 0*(E_om(:,1)-omega_center).^2+0*(E_om(:,1)-omega_center).^3;
[t_tlimit,E_tlimit,omega_tlimit,E_omega_tlimit] = apply_spectral_mask_optimizing(E_om(:,1),E_om(:,2),'1'); % transform-limited case
in.E_om = E_om;
%% Load experimental dScan trace
trace_exp = load('U:\Measurement_Data\OAC_IFROG\2020_10_06_ACER\dscan_combined\data.txt');
wl_exp = load('U:\Measurement_Data\OAC_IFROG\2020_10_06_ACER\dscan_combined\wl.txt');
FS_thickness_exp = load('U:\Measurement_Data\OAC_IFROG\2020_10_06_ACER\dscan_combined\FS_thickness.txt');

in.wl_exp = wl_exp;
in.mat_thickness = FS_thickness_exp-0.2842;
in.trace_exp = trace_exp;
in.material = 'FS';
in.phase_fnctn = 'a(1)*(omega-omega0).^2 + a(2)*(omega-omega0).^3 + a(3)*(omega-omega0).^4';

in.params = [0 0.2;
            -0.05 0.05;
            -0.1 0.1];
in.N_params = [3 3 3];
in.start_params = -1*[-3 +43 +40];
in.start_param_SPM = 0;
in.mask=['exp(1i*(' in.phase_fnctn '))'];
in.flag_plot = 1;
in.flag_progress = 1;
if isempty(in.start_params)
    in.flag_grid = 1;
else
    in.flag_grid = 0;
end
in.seed = [];
%%
out = fit2dScan(in);
%%
retrieved = calc_dscan(out.E_om(:,1),out.E_retrieved,out.material,out.mat_thickness);
load U:\People\S.Peter\projects\IFROG\data_BG39.mat;
data_BG39_ip = interp1(data_BG39(:,1),data_BG39(:,2),retrieved.axis_lambda,'linear');

[t1,Et1] = spectrum2pulse(out.E_om(:,1),out.E_retrieved,'1');

t_cut = 80; % [fs]
T_period = 2*pi/omega_center; % [fs]
points_per_cycle = 10;
Nip = 2^ceil(log2(2*t_cut/T_period*points_per_cycle));
t1_ip = map2colvec(-t_cut:2*t_cut/(Nip-1):t_cut);
if 1
    E1 = flip(Et1,1);
    E1_ip = interp1(t1,conj(E1)/max(abs(conj(E1))),t1_ip,'spline'); % Use 'spline' method for interpolating complex numbers!
else
    E1 = Et1;
    E1_ip = interp1(t1,E1/max(abs(E1)),t1_ip,'spline'); % Use 'spline' method for interpolating complex numbers!
end
pulse_nrg = trapz(t1_ip,abs(E1_ip).^2);
Intensity = abs(E1_ip).^2/pulse_nrg;
dlmwrite('E_vs_t.txt',[map2colvec(t1_ip) map2colvec(real(E1_ip)) map2colvec(imag(E1_ip))],'\t'); % write data to text file
%% Calculate peak intensity

struct_plot.x = out.E_om(:,1);
struct_plot.y1 = abs(out.E_retrieved)/max(abs(out.E_retrieved));
struct_plot.y2 = unwrap(angle(out.E_retrieved));
struct_plot.xlim = [min(out.E_om(:,1)) 3.22];
struct_plot.y1lim = [0 1];
struct_plot.y2lim = [min(unwrap(angle(out.E_retrieved))) max(unwrap(angle(out.E_retrieved)))];
struct_plot.xlabel = '\omega [rad/fs]';
struct_plot.y1_label = 'E(\omega)';
struct_plot.y2_label = '\Phi(\omega)';
struct_plot.title = num2str(out.mask_retrieved);

figure;
subplot(221)
surf(in.wl_exp,in.mat_thickness,in.trace_exp); shading interp; view([0 90]);
xlim([320 500]);
ylim([min(in.mat_thickness) max(in.mat_thickness)]);
xlabel('wavelength [nm]')
ylabel('thickness [mm]')
title('dScan, measured')
colormap jet;
subplot(222)
surf(retrieved.axis_lambda,retrieved.mat_vec,retrieved.trace_lambda.*permute(extend(data_BG39_ip,length(retrieved.mat_vec)),[2 1])); shading interp; view([0 90]);
xlim([320 500]);
ylim([min(in.mat_thickness) max(in.mat_thickness)]);
xlabel('wavelength [nm]')
ylabel('thickness [mm]')
title('dScan, retrieved')

subplot(223)
plotyyP(struct_plot);
setfigP;

subplot(224);
plot(t1_ip,Intensity,'k')
xlabel('time [fs]')
ylabel('|E|^2 (normalized to pulse energy)')
setfigP;
%%
if 0
    save test.mat out;
    dlmwrite('E_vs_omega.txt',[map2colvec(out.E_om(:,1)) map2colvec(real(out.E_retrieved)) map2colvec(imag(out.E_retrieved))],'\t');
    dlmwrite('E_vs_t.txt',[map2colvec(t1_ip) map2colvec(real(E1_ip)) map2colvec(imag(E1_ip))],'\t'); % write data to text file
end
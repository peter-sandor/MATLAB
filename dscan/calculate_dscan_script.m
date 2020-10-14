%% Load Fundamental Spectrum (own measurement)
% filename_sp='U:\Measurement_Data\OAC_IFROG\2019_05_02_ACER\venteon_spectrum_CEP_stabil.txt';
filename_sp='Z:\People\S.Peter\projects\IFROG\Venteon_with_APE_wavescan_noheader.txt';
comma2dot(filename_sp);
venteon_spectrum=load([filename_sp(1:end-4) '_dot.txt']);
wl0=venteon_spectrum(:,1);
sp0=venteon_spectrum(:,2);
% cd U:\Measurement_Data\Nanotriangles_Experiment\20180709_autokorr_tobb_ek_all
sp0=(sp0-mean(sp0(1:50)))/max(sp0-mean(sp0(1:50)));
crop_start=50; crop_end=620;
sp0(sp0<0.01)=0;
% sp0(2910:end)=0;
wl=map2colvec(wl0(crop_start:crop_end));
I_wl=map2colvec(sp0(crop_start:crop_end));
phi_wl = zeros([length(wl) 1]);
%%
E_om = SpectrumConvert([map2colvec(wl) map2colvec(I_wl) map2colvec(phi_wl)]);
[omega_center,sigma_Venteon] = peak_props([E_om(:,1) abs(E_om(:,2)).^2]);
E_om(:,3) = 0*(E_om(:,1)-omega_center).^2+0*(E_om(:,1)-omega_center).^3;
[t_tlimit,E_tlimit,omega_tlimit,E_omega_tlimit] = apply_spectral_mask_optimizing(E_om(:,1),E_om(:,2),'1'); % transform-limited case
%% Calculate phases from material
material = 'FS';
Phi_from_CM = load('Z:\People\S.Peter\projects\IFROG\Layertec_103366_Z0509040_Phi_vs_omega.txt');
L_FS = 1; % [mm]
L_BK7 = 1; % [mm]
L_Air = 1000; % [mm]
nFS = sellmeier(wl*1e-3,'FS');
nBK7 = sellmeier(wl*1e-3,'BK7');
nAir = sellmeier(wl*1e-3,'Air');
% GVD(wl*1e-3,material);
Phi_FS_lambda = 2*pi./(wl*1e-6).*nFS*L_FS;
Phi_BK7_lambda = 2*pi./(wl*1e-6).*nBK7*L_BK7;
Phi_Air_lambda = 2*pi./(wl*1e-6).*nAir.*L_Air;
% figure;plot(wl,Phi_FS_lambda-min(Phi_FS_lambda),'k')
omega0 = map2colvec(2*pi*300./wl);
omega1 = E_om(:,1);
% omega1 = out.data_phasefit(:,1);
Phi_FS = interp1(omega0,Phi_FS_lambda,omega1);
Phi_FS = Phi_FS - (Phi_FS(end)-Phi_FS(1))/(omega1(end)-omega1(1))*omega1;
Phi_FS = Phi_FS - min(Phi_FS);
Phi_BK7 = interp1(omega0,Phi_BK7_lambda,omega1);
Phi_BK7 = Phi_BK7 - (Phi_BK7(end)-Phi_BK7(1))/(omega1(end)-omega1(1))*omega1;
Phi_BK7 = Phi_BK7 - min(Phi_BK7);
Phi_Air = interp1(omega0,Phi_Air_lambda,omega1);
Phi_CM = interp1(Phi_from_CM(:,1),Phi_from_CM(:,2),omega1,'linear','extrap');
Phi_CM = Phi_CM - (Phi_CM(end)-Phi_CM(1))/(omega1(end)-omega1(1))*omega1;
Phi_CM = Phi_CM - min(Phi_CM);
%% Calculate Dscan trace
clear trace_Dscan_lambda trace_Dscan_omega I_SHG_lambda I_SHG_omega;
crystal_length = 0*12; % [um]
crystal_theta = 27.5; % [deg]
mat_length_vec = -2:0.1:2; % [mm]
N_mat = length(mat_length_vec);
phase_mask = 'exp(1i*(0*(omega-omega0).^2 + 20*(omega-omega0).^3 + 0*(omega-omega0).^4))';
% phase_mask = '1';
omega_range = [3.4 6.2];

for ind1 = 1:N_mat
    Phi_total = mat_length_vec(ind1)*Phi_FS;
    [t1,E1,omega,E_omega] = apply_spectral_mask_optimizing(E_om(:,1),E_om(:,2).*exp(1i*Phi_total),phase_mask);
    E_omega_ip = interp1(omega,E_omega,E_om(:,1),'spline');
    [E_omega2,I_lambda] = calc_SHG_spectrum([map2colvec(E_om(:,1)) map2colvec(abs(E_omega_ip)) map2colvec(unwrap(angle(E_omega_ip)))],crystal_length,crystal_theta);
    trace_Dscan_lambda(:,ind1) = I_lambda(:,2);
    trace_Dscan_omega(:,ind1) = E_omega2(:,2);
%     disp([num2str(ind1) '/' num2str(N_mat)]);
end
axis_Dscan_lambda = I_lambda(:,1);
axis_Dscan_omega = E_omega2(:,1);
if 1
    figure;
    surf(axis_Dscan_lambda,mat_length_vec,squeeze(trace_Dscan_lambda(:,:,1)).');
    xlim([300 500]);
    view([0 90]);
    shading interp;
    xlabel('wavelength [nm]');
    ylabel('Path length in FS [mm]');
    colormap jet;
%     figure;surf(lambda_plot_2,mat_length_vec,squeeze(I_SHG_lambda).'); view([0 90]); shading interp; xlabel('wavelength [nm]'); ylabel('glass insertion [mm]');colormap jet;
end
%% Resample time-domain field to obtain arrays of reasonable size
t_cut = 130; % [fs]
T_period = 2*pi/omega_center; % [fs]
points_per_cycle = 10;
Nip = 2^ceil(log2(2*t_cut/T_period*points_per_cycle));
t1 = t_rec2;
t1_ip = map2colvec(-t_cut:2*t_cut/(Nip-1):t_cut);
if 1
    E1 = flip(E_rec2,1);
    E1_ip = interp1(t1,conj(E1)/max(abs(conj(E1))),t1_ip,'spline'); % Use 'spline' method for interpolating complex numbers!
else
    E1 = E_rec2;
    E1_ip = interp1(t1,E1/max(abs(E1)),t1_ip,'spline'); % Use 'spline' method for interpolating complex numbers!
end
% dlmwrite('E_vs_t.txt',[map2colvec(t1_ip) map2colvec(real(E1_ip)) map2colvec(imag(E1_ip))],'\t'); % write data to text file
%% Calculate peak intensity
pulse_nrg = trapz(t1_ip,abs(E1_ip).^2);
Intensity = abs(E1_ip).^2/pulse_nrg;
figure;plot(t1_ip,Intensity,'k')
xlabel('time [fs]')
ylabel('|E|^2 (normalized to pulse energy)')
title(retrieval_path_target)
setfigP;
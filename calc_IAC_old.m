load('U:\Measurement_Data\Hemispherical Autocorrelation Experiment\LaserSpectra\spectrum2018_10_15_17_17.mat');
% cd U:\Measurement_Data\Nanotriangles_Experiment\20180709_autokorr_tobb_ek_all
filename='U:\Measurement_Data\Nanotriangles_Experiment\20180709_autokorr_tobb_ek_all\61p0mm.csv';
sp=(sp-mean(sp(3400:3600)))/max(sp-mean(sp(3400:3600)));
crop_start=1050; crop_end=2950;
sp(sp<0.01)=0;
sp(2910:end)=0;
wl_crop=map2colvec(wl(crop_start:crop_end));
sp_crop=map2colvec(sp(crop_start:crop_end));

A = importdata(filename,',',2);
t=(A.data(:,1)-0.0668)*2.4171/((0.1012-0.028)/15);
om_scorr=FourierAxis(t);
Scorr_meas=(A.data(:,3)-min(A.data(:,3)))/mean(A.data(1:20,3)-min(A.data(:,3)));
Scorr_fft=fft(Scorr_meas);
figure;plot(abs(Scorr_fft),'k')
Scorr_fft_filtered_A=Scorr_fft;
Scorr_fft_filtered_B=Scorr_fft;
Scorr_fft_filtered_C=Scorr_fft;
Scorr_fft_filtered_A(18:584)=0;
Scorr_fft_filtered_B([1:32 59:543 570:600])=0;
Scorr_fft_filtered_C([1:73 102:500 529:600])=0;
Scorr_filtered_A=ifft(Scorr_fft_filtered_A);
Scorr_filtered_B=ifft(Scorr_fft_filtered_B);
Scorr_filtered_C=ifft(Scorr_fft_filtered_C);
figure;plot(t,Scorr_meas,'k')
hold on;plot(t,[Scorr_filtered_A Scorr_filtered_B Scorr_filtered_C])
% figure;plot(fftshift(om_scorr),fftshift(abs(Scorr_fft)))
normfactor0=min(Scorr_filtered_A);
Scorr_meas=Scorr_meas/normfactor0;
t_meas=t;
%% Simulate an ultrashort Gaussian pulse
c=300; % [nm/fs]
tau0=10; % pulse duration, [fs]
omega0=c/800*2*pi;
delta_t=2*pi/omega0/10;
tmax=1/(1/tau0*2*pi/40);
Nt=floor(2*tmax/delta_t);
t=-tmax:2*tmax/(Nt-1):tmax;
field_vs_t=exp(-t.^2/tau0^2).*exp(1i*omega0*t);
field_vs_omega=fft(field_vs_t);
omega=FourierAxis(t);
% figure;plot(t,real(field_vs_t))
% figure;plot(fftshift(omega),fftshift(abs(field_vs_omega)));
[temp1 sigma_om FWHM_om]=calc_stats([map2colvec(fftshift(omega)) map2colvec(abs(fftshift(field_vs_omega)).^2)]);
index_om1=max(vec2ind(omega(1:floor(Nt/2))<=omega0-8*sigma_om));
index_om2=min(vec2ind(omega(1:floor(Nt/2))>=omega0+8*sigma_om));
lambda_to_save=flipdim(c./omega(index_om1:index_om2)*2*pi,2);
intensity_to_save=flipdim(abs(field_vs_omega(index_om1:index_om2)).^2./lambda_to_save.^2,2)/max(abs(field_vs_omega(index_om1:index_om2)).^2./lambda_to_save.^2);
wl=lambda_to_save;
sp=intensity_to_save;
dlmwrite('simulated_spectrum.txt',[map2colvec(lambda_to_save) map2colvec(intensity_to_save)],'\t');
%
tau=t;
t_meas=t;
Et1=field_vs_t;
clear Scorr;
for ind1=1:Nt
    index_start=min(vec2ind(t>=max(t(1),t(1)+tau(ind1))));
    index_end=max(vec2ind(t<=min(t(end),t(end)+tau(ind1))));
    Et2=interp1(t,Et1,t(index_start:index_end)-tau(ind1),'linear');
    Scorr(ind1)=sum(abs((Et1(index_start:index_end)+Et2).^2).^2);
    Scorr2a(ind1)=sum(abs(Et1(index_start:index_end)).^4)+sum(abs(Et2).^4); % constant term
    Scorr2b(ind1)=4*sum(abs(Et1(index_start:index_end)).^2.*abs(Et2).^2); % intensity autocorrelation
    Scorr2c(ind1)=4*sum((abs(Et1(index_start:index_end)).^2+abs(Et2).^2).*real(Et1(index_start:index_end).*conj(Et2))); % interferometric term at the fundamental frequency
    Scorr2d(ind1)=2*sum(real(Et1(index_start:index_end).^2.*conj(Et2).^2)); % interferometric term at the second harmonic
end
Scorr2=Scorr2a+Scorr2b+Scorr2c+Scorr2d;
index_norm=min(vec2ind(tau>=(tau(1)+1.5*tau0)));
% normfactor=Scorr(index_norm);
normfactor=Scorr2a(round(Nt/2));
Scorr=Scorr/normfactor;
Scorr2=Scorr2/normfactor;
Scorr_meas=Scorr;
figure;hold on;
plot(tau,Scorr,'k')
plot(tau,Scorr2,'r')
plot(tau,[Scorr2a.' Scorr2b.' Scorr2c.' Scorr2d.']/normfactor);
wl_crop=wl;
sp_crop=sp;
%%
sp_omega=SpectrumConvert([map2colvec(wl_crop) map2colvec(sp_crop)]);
% [omega0 sigma_om FWHM]=calc_stats([map2colvec(sp_omega(:,1)) abs(sp_omega(:,2)).^2]);

in.t_meas=t_meas;
in.Scorr_meas=Scorr_meas;
in.mask0='1';
in.sp_omega=sp_omega;
[taxis td_data om sp_om]=simulate_pulse_shaper(sp_omega,in.mask0);
%%
out = fit2IAC(in);
record_sorted=sortrows(out.record,3);
record_sorted(1,:)

Eomega=fft(Et1);
figure;
subplot(221);plot(wl_crop,sp_crop,'k');ylim([0 max(sp_crop)]);xlabel('wavelength [nm]')
subplot(222);plot(sp_omega(:,1),sp_omega(:,2),'k');xlabel('\omega [2*pi/fs]')
subplot(223);plot(t,Et1,'k');xlabel('delay [fs]');ylabel('E (normalized)')
subplot(224);
hold on
plot(t,Scorr_meas,'k')
plot(t,Scorr,'r');
ylim([min(Scorr/mean(Scorr(1:100))) max(Scorr/mean(Scorr(1:100)))])
xlabel('delay [fs]')
ylabel('AC signal (normalized)')
load wavelengths.txt;
load delays.txt;
load FROG.txt;

N=256;
freq0=flipdim(299.792./wavelengths,2);
T0=5.0333;
delays_fs0=(delays-T0)*1e6/299.792*2;
freq_center=freq0(round(length(freq0)/2));
freq_step=min([(freq_center-min(freq0))/(N/2-1) (max(freq0)-freq_center)/(N/2-1)]);
freq_min=freq_center-(N/2-1)*freq_step;
freq_max=freq_center+N/2*freq_step;

freq = (freq_min:freq_step:freq_max);
freq = freq - freq(round(N/2)) + freq_center;

delay_step=(max(delays_fs0)-min(delays_fs0))/(N-1);
delays_fs=min(delays_fs0):delay_step:max(delays_fs0);
FROG_freq=flipdim(FROG.',1);

[DELAY0,FREQ0] = meshgrid(delays_fs0,freq0);
[DELAY,FREQ] = meshgrid(delays_fs,freq);
FROG_ip = interp2(DELAY0,FREQ0,FROG_freq,DELAY,FREQ);
FROG_ip=FROG_ip/max(max(FROG_ip));

% figure;surf(delays_fs0,freq0,FROG_freq)
% shading interp;
% view([0 90]);
% 
figure;surf(delays_fs,freq,FROG_ip)
shading interp;
view([0 90]);
% seed=exp();
units{1}='fs';
units{2}='PHz';
%% calculate seed from measured spectrum
spectrum=spproc2('I:\data2017\2017_02_13_PLSB1\amp_spectrum.txt');
sp_omega=SpectrumConvert(spectrum);
seed0=simulate_pulse_shaper([sp_omega(:,1) sp_omega(:,2)],'1'); % sp_omega and seed contains the electric fields as a function of omega and time, respectively.
seed=map2colvec(interp1(seed0(:,1),seed0(:,2),delays_fs));
seed(isnan(seed))=0;
% [fnlimg fnldt] = prepFROG(delay_step,freq_step,N,0,1,0,0);
%%
[Pt, Fr, G, iter] = svdFROG(FROG_ip, seed, 1e-2, 5000,'[1,0,1,0]',0, delay_step, units);
figure;imagescP(delays_fs,freq,Fr)
figure;hold on;
plot(delays_fs,abs(seed),'k')
plot(delays_fs,abs(Pt)/max(abs(Pt)),'r')
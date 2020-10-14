N=2048;
omega_max=2.8;
omega_min=-omega_max;
%%
N=256;
omega_max=2.8;
omega_min=2.0;
%% define spectral amplitude and phase a.f.o. omega
omega0=2*pi*0.3846;
% omega0=2*pi*0.05;
sigma_omega=0.12;
omega_step=(omega_max-omega_min)/(N-1);
omega=omega_min:omega_step:omega_max;
amp_sp=exp(-2.355*(omega-omega0).^2/sigma_omega^2);
% M='exp(1i*0*(omega-omega0).^2)';
M='exp(1i*1*(1000)*(omega-omega0).^3)';
phase_sp=unwrap(angle(eval(M)));
field1=amp_sp.*exp(1i*phase_sp);
% figure;plot(omega,amp_sp,'k');
% hold on;plot(omega,unwrap(angle(field1)),'b')
% generate FROG trace
[data_td] = simulate_pulse_shaper([map2colvec(omega) map2colvec(amp_sp)],M);
timeax=data_td(:,1);
field_td=data_td(:,2);
[omegaax,tauax,signal_SHG,signal_SD] = FROG_calc(timeax,field_td);

FROG_SD=fftshift(signal_SD,1);
FROG_SHG=fftshift(signal_SHG,1);
omegaaxx=fftshift(omegaax);
% figure;imagescP(tauax,omegaaxx,FROG_SD)
% figure;imagescP(tauax,omegaaxx,FROG_SHG)
% convert FROG trace from omega to lambda-dependence

[temp1,temp2,temp3,temp,crop_SD_delay]=calc_stats([tauax map2colvec(squeeze(sum(FROG_SD,1)))/max(squeeze(sum(FROG_SD,1)))]);
[temp1,temp2,temp3,temp,crop_SD_omega]=calc_stats([omegaaxx map2colvec(squeeze(sum(FROG_SD,2)))/max(squeeze(sum(FROG_SD,2)))]);
[temp1,temp2,temp3,temp,crop_SHG_delay]=calc_stats([tauax map2colvec(squeeze(sum(FROG_SHG,1)))/max(squeeze(sum(FROG_SHG,1)))]);
[temp1,temp2,temp3,temp,crop_SHG_omega]=calc_stats([omegaaxx map2colvec(squeeze(sum(FROG_SHG,2)))/max(squeeze(sum(FROG_SHG,2)))]);

crop_SD_omega(1)=crop_SD_omega(1)+0*20;
FROG_SHG_cropped=FROG_SHG(crop_SHG_omega(1):crop_SHG_omega(2),crop_SHG_delay(1):crop_SHG_delay(2));
FROG_SD_cropped=FROG_SD(crop_SD_omega(1):crop_SD_omega(2),crop_SD_delay(1):crop_SD_delay(2));
tau_SHG=tauax(crop_SHG_delay(1):crop_SHG_delay(2));
tau_SD=tauax(crop_SD_delay(1):crop_SD_delay(2));
omega_SHG=omegaaxx(crop_SHG_omega(1):crop_SHG_omega(2));
omega_SD=omegaaxx(crop_SD_omega(1):crop_SD_omega(2));
N2=256;
c=299.792;
%% omega --> lambda switch for SHG trace (equidistant lambda)
lambda_SHG=2*pi*c./omega_SHG;
lambda_SHG_min=min(lambda_SHG);
lambda_SHG_max=max(lambda_SHG);
lambda_SHG_span=lambda_SHG_max-lambda_SHG_min;
lambda_SHG_step=lambda_SHG_span/(N2-1);
lambda_SHG_ip=lambda_SHG_min:lambda_SHG_step:lambda_SHG_max;
tau_SHG_ip=min(tau_SHG):(max(tau_SHG)-min(tau_SHG))/(N2-1):max(tau_SHG);
FROG_SHG_lambda=FROG_SHG_cropped./extend(lambda_SHG,length(tau_SHG)).^2;
[TAU_SHG,LAM_SHG]=meshgrid(tau_SHG,lambda_SHG);
[TAU_SHG_ip,LAM_SHG_ip]=meshgrid(tau_SHG_ip,lambda_SHG_ip);
FROG_SHG_ip = interp2(TAU_SHG,LAM_SHG,FROG_SHG_lambda,TAU_SHG_ip,LAM_SHG_ip);
delay_SHG_mm=tau_SHG_ip*299.792/1e6/2;
dlmwrite('FROG.txt',FROG_SHG_ip/max(max(FROG_SHG_ip)),'delimiter','\t','precision','%.3E');
dlmwrite('delays.txt',delay_SHG_mm,'\t');
dlmwrite('wavelengths.txt',lambda_SHG_ip,'\t');
FROG2ifa('.\FROG.txt');
figure;imagescP(tau_SHG_ip,lambda_SHG_ip,FROG_SHG_ip)
%% omega --> lambda switch for SD trace (equidistant lambda)
lambda_SD=2*pi*c./omega_SD;
lambda_SD_min=min(lambda_SD);
lambda_SD_max=max(lambda_SD);
lambda_SD_span=lambda_SD_max-lambda_SD_min;
lambda_SD_step=lambda_SD_span/(N2-1);
lambda_SD_ip=lambda_SD_min:lambda_SD_step:lambda_SD_max;
tau_SD_ip=min(tau_SD):(max(tau_SD)-min(tau_SD))/(N2-1):max(tau_SD);
FROG_SD_lambda=FROG_SD_cropped./extend(lambda_SD,length(tau_SD)).^2;
[TAU_SD,LAM_SD]=meshgrid(tau_SD,lambda_SD);
[TAU_SD_ip,LAM_SD_ip]=meshgrid(tau_SD_ip,lambda_SD_ip);
FROG_SD_ip = interp2(TAU_SD,LAM_SD,FROG_SD_lambda,TAU_SD_ip,LAM_SD_ip);
delay_SD_mm=tau_SD_ip*299.792/1e6/2;
dlmwrite('FROG.txt',FROG_SD_ip/max(max(FROG_SD_ip)),'delimiter','\t','precision','%.3E');
dlmwrite('delays.txt',delay_SD_mm,'\t');
dlmwrite('wavelengths.txt',lambda_SD_ip,'\t');
FROG2ifa('.\FROG.txt');
figure;imagescP(tau_SD_ip,lambda_SD_ip,FROG_SD_ip)
%% omega --> lambda switch for SD trace (nonequidistant lambda, equidist. omega)
omega_SD_ip=min(omega_SD):(max(omega_SD)-min(omega_SD))/(N2-1):max(omega_SD);
tau_SD_ip=min(tau_SD):(max(tau_SD)-min(tau_SD))/(N2-1):max(tau_SD);
[TAU_SD,OMEGA_SD]=meshgrid(tau_SD,omega_SD);
[TAU_SD_ip,OMEGA_SD_ip]=meshgrid(tau_SD_ip,omega_SD_ip);
FROG_SD_ip = interp2(TAU_SD,OMEGA_SD,FROG_SD_cropped,TAU_SD_ip,OMEGA_SD_ip);
delay_SD_mm=tau_SD_ip*299.792/1e6/2;
lambda_SD_nonequiv=2*pi*c./flipdim(map2colvec(omega_SD_ip),1);
dlmwrite('FROG.txt',flipdim(FROG_SD_ip/max(max(FROG_SD_ip)),1),'delimiter','\t','precision','%.3E');
dlmwrite('delays.txt',delay_SD_mm,'\t');
dlmwrite('wavelengths.txt',lambda_SD_nonequiv,'\t');
FROG2ifa('.\FROG.txt');
% figure;surf(tau_SD_ip,lambda_SD_nonequiv,flipdim(FROG_SD_ip,1));shading flat;view([0 90])
%% convert to spectrum a.f.o lambda
[temp1,temp2,temp3,temp,crop_omega]=calc_stats([map2colvec(omega) map2colvec(amp_sp)/max(amp_sp)]);
omega_crop=omega(crop_omega(1):crop_omega(2));
lambda=2*pi*c./omega_crop;
lambda_min=min(lambda);
lambda_max=max(lambda);
lambda_span=lambda_max-lambda_min;
lambda_step=lambda_span/(N2-1);
lambda_ip=lambda_min:lambda_step:lambda_max;
spectrum=map2rowvec(amp_sp(crop_omega(1):crop_omega(2))).^2./lambda.^2; % apply Jacobian: f(lambda)=g(omega)*omega^2=g(2*pi*c/lambda)/lambda^2
spectrum=spectrum/max(spectrum);
spectrum_ip=interp1(lambda,spectrum,lambda_ip,'linear');

figure;plot(lambda_ip,spectrum_ip,'k')
dlmwrite('spectrum.txt',[map2colvec(lambda_ip) map2colvec(spectrum_ip)],'\t');
plotreconst2('spectrum.txt');fitreconst;
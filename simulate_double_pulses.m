% calib=load('L:\2012 05 22 E\uvAOMcalibration');
% [fs]
Nk=500;
omegaCenter=2*pi*1.155;
omega_span=0.05;
omega_min=omegaCenter-omega_span/2*pi;
omega_max=omegaCenter+omega_span/2*pi;
omega_step=(omega_max-omega_min)/(Nk-1);
omega=(omega_min:omega_step:omega_max); % dimension is thought to be in [pixels], but it could as well be in [ns]
t_step=2*pi/(omega_max-omega_min);
t=-Nk/2*t_step+t_step:t_step:+Nk/2*t_step;

% omegaR=2*pi*1.155e6; % [rad/ns]
omegaR=0; % [rad/ns]
sigma=omega_span/2; 
% tR=(omegaR-calib(2))/calib(1);
tR=0;
N=4;
figure('color','w')
for ind1=1:N
    ppdelay=ind1*125; % [fs]
    phase=(omega-omegaR)*ppdelay;
    % phase2=sin(calib(1)*(tAC-tR)*ppdelay);
    sp_opt=1/2*exp(-(omega-omegaCenter).^2/sigma^2).*(exp(i*phase)+1);
    % sp_opt=1/2*exp(-(tAC-tCenter).^2/sigma^2).*exp(i*phase2);
    % phase=unwrap(angle(sp_opt));
    pulses=ifftshift(ifft([sp_opt, zeros([1, 3*Nk])]));
    hndl(ind1,1)=subplot(N,2,(ind1-1)*2+1);
    plot(abs(pulses).^2,'k','linewidth',3);
%     xlabel('delay [fs]')
    % hold on;plot(tAC,real(pulses),'r')
    % plot(tAC,imag(pulses),'b')
    % hold off;
    xlim([920 1020])
    hndl(ind1,2)=subplot(N,2,(ind1-1)*2+2);
    plot(omega,abs(sp_opt).^2,'k','linewidth',3);
%     xlabel('frequency [rad/fs]')
    xlim([7.2 7.3])
    % delay in [fs] from graph: (delay read off in pixels)/(3360*calib(1))*2*pi
end
set(hndl,'Visible','off')
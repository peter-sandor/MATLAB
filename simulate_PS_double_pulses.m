calib=load('L:\2012 05 22 E\uvAOMcalibration');
sigma=300; % [pixels]
tAC=(1:3360); % dimension is thought to be in [pixels], but it could as well be in [ns]
%%
tCenter=1680;
% omegaR=2*pi*1.155e6; % [rad/ns]
omegaR=2*pi*0; % [rad/ns]
tR=(omegaR-calib(2))/calib(1);

figure('color','w')
for ind1=1:3
    ppdelay=ind1*6.9524e-4; % [ns]
    phase=calib(1)*(tAC-tR)*ppdelay;
    % phase2=sin(calib(1)*(tAC-tR)*ppdelay);
    sp_opt=1/2*exp(-(tAC-tCenter).^2/sigma^2).*(exp(i*phase)+1);
    % sp_opt=1/2*exp(-(tAC-tCenter).^2/sigma^2).*exp(i*phase2);
    % phase=unwrap(angle(sp_opt));
    pulses=ifftshift(ifft(sp_opt));
    hndl(ind1,1)=subplot(3,2,(ind1-1)*2+1);
    plot(abs(pulses));
    xlabel('delay [arb. units]')
    % hold on;plot(tAC,real(pulses),'r')
    % plot(tAC,imag(pulses),'b')
    % hold off;
    xlim([1610 1700])
    hndl(ind1,2)=subplot(3,2,(ind1-1)*2+2);
    plot(abs(sp_opt).^2,'k');
    xlabel('frequency [arb. units]')
    xlim([1000 2300])
    % delay in [ns] from graph: (delay read off in pixels)/(3360*calib(1))*2*pi
end
set(hndl,'Visible','off')
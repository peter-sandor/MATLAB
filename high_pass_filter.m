R=1e4;
C=2e-9;
omega_3dB=1/R/C;
nu_3dB=omega_3dB/2/pi;
Nomega=100;
omega_max=omega_3dB*5;
omega=(0:omega_max/(Nomega-1):omega_max);
Z=R./sqrt(R^2+(1./omega.^2/C^2)).*exp(-1i*atan(1./omega/R/C));
figure;
subplot(211)
plot(omega/2/pi,abs(Z),'k')
xlabel('\nu [Hz]')
ylabel('|Z|')
title(['\nu_{3dB} = ' num2str(roundP(nu_3dB,0)) ' Hz'])
subplot(212)
plot(omega/2/pi,unwrap(angle(Z))/pi,'b')
xlabel('\nu [Hz]')
ylabel('phase [\pi]')
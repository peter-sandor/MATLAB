uiopen('I:\users\ps5at\DAMOP_2017\figures\OCS_3panel_data_only.fig',1)
[delay,yield]=getdata;
%%
Ndelay=128;
Next=0*Ndelay;
delay_min=70;
delay_max=max(delay);
delay_step=(delay_max-delay_min)/(Ndelay-1);
delay_ip=delay_min:delay_step:delay_max;
yield_ip=interp1(delay,yield,delay_ip);
delay_ext=[(delay_ip(1)-Next*delay_step):delay_step:(delay_ip(1)-delay_step) delay_ip];
yield_ext=[zeros([1 Next]) yield_ip];
omega_ax=FourierAxis(delay_ext);
yield_fft=fft(yield_ext-mean(yield_ext));
figure;plot(fftshift(omega_ax)/2/pi,fftshift(abs(yield_fft)),'k')

%%
% rotational constants
B=0.203; % [1/cm]
D=3.46e-8; % [1/cm]
T=10;

maxJ=20;
J=0:maxJ;
h=6.626e-34; % [J*s]
c=3e8; % [m/s]
E=h*c*100*(B*J.*(J+1)+D*J.^2.*(J+1).^2);
p=MB_distr(E,T);
coeffsJ=p/sum(p);

% E_sp=fftshift(omega_ax)*h/2/pi*1e12/1.6e-19; % [eV]
J_omega=2*pi*3e-2*(B*J.*(J+1)+D*J.^2.*(J+1).^2); % [2*pi*THz]
J_bin_edges(1)=0;
for ind1=2:length(J_omega)-1
    J_bin_edges(ind1)=mean([J_omega(ind1) J_omega(ind1+1)]);
end
J_omega_ip=sort([omega_ax(1:Ndelay) J_omega J_bin_edges]);
sp_ip=interp1(fftshift(omega_ax),fftshift(abs(yield_fft)),J_omega);
sp_ip2=interp1(fftshift(omega_ax),fftshift(abs(yield_fft)),J_omega_ip);
sp_integr(1)=J_bin_edges(1)*interp1(fftshift(omega_ax),fftshift(abs(yield_fft)),J_bin_edges(1));
for ind1=1:length(J_bin_edges)-1
    domain_x(:,ind1)=map2colvec(J_bin_edges(ind1):(J_bin_edges(ind1+1)-J_bin_edges(ind1))/(9):J_bin_edges(ind1+1));
    domain_y(:,ind1)=map2colvec(interp1(fftshift(omega_ax),fftshift(abs(yield_fft)),domain_x(:,ind1)));
    sp_integr(ind1+1)=trapz(domain_x(:,ind1),domain_y(:,ind1));
end
sp_integr(1)=sp_integr(1)*2;
figure;
subplot(2,2,[1 2])
plot(delay,yield,'k.--')
xlabel('delay [ps]')
ylabel('yield')
subplot(223)
plot(fftshift(omega_ax)/2/pi*1e3,fftshift(abs(yield_fft)),'k')
hold on;
plot(J_omega/2/pi*1e3,sp_ip,'r*')
% plot(J_omega_ip/2/pi*1e3,zeros([length(J_omega_ip) 1]),'b*')
% plot(J_omega(1:end-1)/2/pi*1e3,sp_integr,'bo')
xlim([0 max(fftshift(omega_ax)/2/pi*1e3)])
xlabel('\nu [GHz]')
subplot(224)
bar(J,coeffsJ)
hold on;
bar(J(1:length(sp_integr)),sp_integr/sum(sp_integr),'r')
xlabel('J')
xlim([0 maxJ])
title(['T = ' num2str(T) ' K'])
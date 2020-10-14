%Everything in Atomic units unless otherwise specified
clear
T=1.0*2.7/0.024; % optical period
tau=10*T; % temporal extent of field [24 as]
dt=T/20;
tmin=0;
tmax=3*tau;
t0=0;%tmax/2;
t=tmin:dt:tmax;
% for Ti:Sapphire frequency which has a 2.7 fs period
omega=2*pi/T;
sigma=1*20/0.53e-4; % spatial extent of field in transverse direction [Bohr]
phi=0*pi/2;
init_cond=[100,0]; % [x dx/dt]
% Calculation for pulse intensity of 1e13 W/cm^2, which corresponds to a field of 2.74e9 V/m, which is  
% 0.0053 atomic units.  Thus, the vector potential is 0.0053/omega
field_amp=+1e-1*0.0053;
envelope=-field_amp/omega*exp(-(t-t0).^2/tau^2);
A=envelope.*cos(omega*t+phi);
%Note that in Atomic units, Up=A^2/4 (factor of two for cycle average and
%take energy, which goes like A^2.
up=envelope.^2/4;
% func_e_in_field = @(t,x) [x(2); field_amp*exp(-x(1)^2/sigma^2)*exp(-(t-t0)^2/tau^2)*cos(omega*t+phi)];
func_e_in_field = @(t,x) [x(2); field_amp*x(1)*heaviside(x(1))*cos(omega*t+phi)];
[tout, yout] = ode45(func_e_in_field,[tmin tmax],init_cond);
yip=interp1(tout,yout,t);

temp=zeros([2 length(t)]);
for ind1=1:length(t)
    temp(:,ind1)=func_e_in_field(t(ind1),init_cond);
end
field_time=temp(2,:);

temp=zeros([2 length(t)]);
for ind1=1:length(t)
    temp(:,ind1)=func_e_in_field(t(ind1),[yip(ind1,1) 0]);
end
field_actual=temp(2,:);

% field_time2=field_amp.*exp(-(map2colvec(t)-t0).^2/tau^2).*cos(omega*map2colvec(t)+phi);
% field_actual2=field_amp*exp(-yip(:,1).^2/sigma^2).*exp(-(map2colvec(t)-t0).^2/tau^2).*cos(omega*map2colvec(t)+phi);
for ind1=1:length(t)
    E_timeint=sum(exp(-(map2colvec(t(1:ind1))-t0).^2/tau^2).*cos(omega*map2colvec(t(1:ind1))+phi))*dt;
    A_actual(ind1)=field_amp*exp(-yip(ind1,1).^2/sigma^2).*E_timeint;
end
up2=A_actual.^2/4;

Ekin=1/2*(yip(:,2)).^2;
steps_to_avg=ceil(2*pi/omega/dt);
for ind1=1:length(Ekin)-steps_to_avg
    temp=circshift(Ekin,-[ind1-1 0]);
    Ekin_cycavg(ind1)=sum(temp(1:steps_to_avg))/steps_to_avg;
end
% xmax=max(max(abs(yout(:,1))),2*sigma);
xmax=max(abs(yip(:,1)));
Nx=100;
xvec=-xmax:2*xmax/(Nx-1):xmax;
field_spatial_profile=field_amp*exp(-xvec.^2/sigma^2);
max(xvec)/sigma;
%%
figure(1)
ax(1)=subplot(221);
plot(t,yip(:,1),'k')
hold on;
plot(t,init_cond(1)-1/pi^2*4*field_amp^2/omega^2*init_cond(1)/2.*map2colvec(t).^2,'r')
hold off;
xlabel('time [atomic units]')
ylabel('x [Bohr]')
xlim([tmin tmax]);
ylim([0 +xmax]);
% ylim([-xmax +xmax]);
% ax(2)=subplot(222);
% plot(field_spatial_profile,xvec,'r')
% ylabel('x [Bohr]')
% xlabel('field strength');
% ylim([-xmax +xmax]);
% xlim([0 max(field_spatial_profile)])
ax(3)=subplot(223);
plot(t,field_time,'r')
hold on
plot(t,field_actual,'b--')
hold off
ylabel('field strength');
legend('E_{pulse}','E @ electron location')
xlim([tmin tmax]);
ax(4)=subplot(224);
hold off;
plot(t,Ekin-1/2*init_cond(2)^2,'k')
hold on
plot(t(1:end-steps_to_avg),Ekin_cycavg-1/2*init_cond(2)^2,'r')
plot(t,up2,'g')
plot(t,up,'m')
xlabel('time [atomic units]')
ylabel('Energy [at. units]')
legend('E_{kin}-E_0','E_{kin} cycle-averaged -E_0','A_{@electron}^2/4','A_{envelope}^2/4')
xlim([tmin tmax]);
% linkaxes(ax([1 3 4]),'x')
% linkaxes(ax([1 2]),'y')

%%
% n=0;
% x=0:xmax/(Nx-1):xmax;
% % F0=F0_max*2/sigma^2*x.*exp(-x.^2/sigma^2);
% F0=field_amp*x;
% index=Nx;
% for ind1=index
%     q=2*F0(ind1)/(omega^2);
%     [q,n]= meshgrid(q,n);
%     [ce,lc]= cen(omega*t,q,n);
%     amp(ind1)=max(ce)-min(ce);
% end
% % [se,ls]= sen(x,q,n);
% 
% figure(2)
% subplot(121)
% plot(t,ce/mean(ce)*x(index))
% xlim([0 5*T])
% subplot(122)
% plot(x,F0,'k')
% % hold on;
% % plot(x,amp,'r') 
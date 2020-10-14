function moADK
% define pulse shape in atomic units:
% [t]=24.1888 attosec
% [E]=5.1422e11 V/m
% [I]=6.43641e15 W/cm^2
t_min=0;
t_max=7000;
t_step=10;
taxis=map2colvec(t_min:t_step:t_max);
Nt=length(taxis);
Topt=107; % 2.6 fs period --> lambda=800 nm
% 1e13 - 1e15 W/cm^2 intensity --> 8.68e+08 - 8.68e+9 V/m field
field_exp=7:0.05:11.8;
F0=10.^field_exp/5.14e11; % field in atomic uits
I0=1/2*(F0*5.14e11).^2*3e8*8.85e-12./1e4; % W/cm2
sigma1=1529; % FWHM=37 fs pulse, sigma=FWHM/sqrt(2*ln(2))
sigma2=1240;
delay1=t_max/2;
delay2=delay1+2150+250;	% 52+1.3 fs delay
phase1=0;
phase2=0:pi/10:2*pi;
omega=2*pi/Topt;
field=map2colvec(exp(-(taxis-delay1).^2/sigma1^2).*cos(omega*(taxis-delay1)+phase1))*map2rowvec(F0);

intensity=137/8/pi*(field).^2;
T=taxis;
F=abs(field);

%charges for the two 'nuclei' in our 'diatomic molecule'
Za=1;
Z=0;

m=0;
l=0;
l1=1;
l2=2;
l3=3;
l4=4;
am=abs(m);
molecule='N2+';

switch molecule
    case 'H2+'
        IP=0/27.211;
        Cl0=4.37;
        Cl1=0;
        Cl2=0.05;
        Cl3=0;
        Cl4=0;
    case 'N2+'
        IP=15.58/27.211;
        Cl0=2.02;
        Cl1=0;
        Cl2=0.78;
        Cl3=0;
        Cl4=0.04;
    case 'SO+'
        IP=0/27.211;
        Cl0=0;
        Cl1=0.41;
        Cl2=-0.31;
        Cl3=0.01;
        Cl4=0;
    case 'NO(pi)'
        IP=0/27.211;
        Cl0=0;
        Cl1=0.22;
        Cl2=1.15;
        Cl3=0.01;
        Cl4=0.0;
%       R0=0;
end

kappa=sqrt(2*IP);
%factorial term
q=(-1)^m*sqrt(((2*l+1)*factorial(l+am))/(2*factorial(l-am)));
q1=(-1)^m*sqrt(((2*l1+1)*factorial(l1+am))/(2*factorial(l1-am)));
q2=(-1)^m*sqrt(((2*l2+1)*factorial(l2+am))/(2*factorial(l2-am)));
q3=(-1)^m*sqrt(((2*l3+1)*factorial(l3+am))/(2*factorial(l3-am)));
q4=(-1)^m*sqrt(((2*l4+1)*factorial(l4+am))/(2*factorial(l4-am)));
%B^2term term
B2=((Cl0*q)+(Cl1*q1)+(Cl2*q2)+(Cl3*q3)+(Cl4*q4))^2;
%B2=((Cl.*q)).^2;
%Now the tunnelling ionization rate in a static field (= wstat)
wstat_DC=(B2/((2^am)*factorial(am))).*(1/(kappa^((2*Z)/kappa-1))).*((2*kappa^3)./F).^((2*Z)/kappa-am-1).*exp((-2*kappa^3)./(3.*F));
wstat_AC=sqrt(3*permute(extend(F0,Nt),[2 1])/pi/kappa).*wstat_DC;
wprob=1-exp(-1*sum(abs(wstat_AC),1)*t_step);

PREC_ABS = 5E-4; %Absolute Precision
PREC_REL = 5E-4; %Relative Precision
options = odeset('RelTol',PREC_ABS,'AbsTol',PREC_REL); %Options for integrator
for ind1=1:size(wstat_AC,2)
    [t_out,pops] = ode45(@RateEq,[t_min t_max],[t_min;1;0],options);
    pops_ip(:,1,ind1)=interp1(t_out,pops(:,1),taxis);
    pops_ip(:,2,ind1)=interp1(t_out,pops(:,2),taxis);
    pops_ip(:,3,ind1)=interp1(t_out,pops(:,3),taxis);
    yield_total(ind1)=squeeze(pops_ip(end,end,ind1));
    [vol1 vol2]=VolAvg(I0(1:ind1)/I0(ind1));
    yield_vavg(ind1)=sum(vol1.*yield_total(1:ind1));
end
% yield_total=squeeze(pops_ip(end,end,:));

figure
subplot(122)
loglog(I0,yield_total,'ko-')
hold on;
loglog(I0,wprob,'bo-')
% loglog(I0,yield_vavg,'mo-')
setfigP
%axis([0.05 1 1e6 1e18])
subplot(121)
hold on;
plot(I0,yield_total,'ko-')
plot(I0,wprob,'bo-')
% plot(I0,yield_vavg,'mo-')
setfigP

function pops_dot = RateEq(t,pops)
    wstat_ip=interp1(taxis,wstat_AC(:,ind1),t);
    pops_dot(1)=1; % this variable is time
    pops_dot(2)=-wstat_ip*pops(2)*pops_dot(1); % population of bound state
    pops_dot(3)=wstat_ip*pops(2)*pops_dot(1); % population of continuum
    pops_dot=map2colvec(pops_dot);
end
end
% Calculate properties of molecular beams
% Based on Moore, Davis, Coplan: Building Scientific Apparatus (4th ed.)
% Chapter 3 (subsection 3.5.6)

T_ratio=273; % [K]
P_Torr=1e-3;
P_SI=P_Torr/750*1e5;
kB=1.38e-23; % [J/K]
density=P_SI/kB/T_ratio; % [1/m^3]
mass=2*1e-27; % [kg]
mdiam=106e-12; % [m]

v_bar=sqrt(8*kB*T_ratio/pi/mass);
Z=sqrt(2)*pi*density*mdiam^2*v_bar;
mean_free_path=1/(sqrt(2)*pi*density*mdiam^2);

%% expansion from a nozzle

length=20-3;
diameter=150e-6;
M=mass*6e26; % [g/mole]
P_back_Torr=750;
P_back_SI=P_back_Torr/750*1e5;

% conversion to reduced pressure
P_R=P_back_SI*(length/1e-2)*(mdiam/1e-10)^2;
% calculate reduced beam width, flux and throughput
H_R=2.48e2*sqrt(P_R);
I_R=1.69e20*sqrt(P_R);
q_R=2.16e21*P_R;

H=H_R/(length/1e-2)*(diameter/1e-2);
I=sqrt(T_ratio/(295*M))*(diameter/1e-2)^2*I_R/((length/1e-2)*(mdiam/1e-10)^2);
q=sqrt(T_ratio/(295*M))*(diameter/1e-2)^3*q_R/((length/1e-2)^2*(mdiam/1e-10)^2);

%%
% the following is from Giacinto Scoles: Atomic and Molecular Beam Methods
% V1 (1988) Sections 2.1.2 and 2.2.3


R=8.314; % [J/K/mole]
P0=10*1e5; % stagnant (backing) pressure [Pa]
Pb=1e-8/760*1e5; % background pressure [Pa]
gamma=7/5;
W=60e-27*6e23; % molar weight
d=150e-6; % nozzle diameter
z=5:5:100; % axial distance in nozzle diameters (d)
z0=0.40;
T0=300; % [K]
A=3.65;
phi=1.662;
theta=0:pi/2/49:pi/2;

zM=0.67*sqrt(P0/Pb); % Mach disk position in nozzle diameters
M=A*(z-z0).^(gamma-1)-1/2*(gamma+1)/(gamma-1)/A./(z-z0).^(gamma-1); % Mach number as a function of distance from nozzle
T_ratio=1./(1+(gamma-1)/2.*M.^2); % relative temperature (compared to stagnant temperature T0)
v_sound=sqrt(gamma*R*T_ratio*T0/W);
v=M.*sqrt(gamma*R*T0/W)./sqrt(1+(gamma-1)/2*M.^2);
P_ratio=T_ratio.^(gamma/(gamma-1));
rho_ax_ratio=T_ratio.^(1/(gamma-1)); % density on axis
rho_offax_ratio=cos(pi*theta/2/phi).^2; % density dependence on azimuthal angle

figure;hold on;
plot(z,M,'k.-');
plot(z,T_ratio,'b.-');
plot(z,v,'r.-');
plot(z,P_ratio,'g.-');
plot(z,rho_ax_ratio,'m.-');
constants_si

FE = 11.8; % Field Enhancement
WorkFunction = 5.3; % [eV]

FWHM = 1e-3; % [m]
lambda = 800e-9; % [m]
% intensity FWHM  = sqrt(2*log(2))*w, where Intensity = I0*exp(-2*x^2/w^2)
w = 0.85*FWHM;
f = 50e-3; % [m]
% w0=lambda*f/(pi*w);
w0x = 5e-6; % [m]
w0y = 5e-6; % [m]
focal_area = pi*w0x*w0y; % [m^2]
beam_power = 1e-3; % [W]
intensity=beam_power/focal_area; % [W/m^2]
reprate=8e7; % [Hz]
pulse_duration = 7e-15; % [s]
pulse_energy = beam_power/reprate; % [J]
% pulse_energy = 1e-9;
peak_intensity = pulse_energy/pulse_duration/focal_area; % [W/m^2]
% peak_intensity = 1e15; % [W/m^2]
field = sqrt(2*peak_intensity/C_SI.eps0/C_SI.c); % [V/m]
omega = 2*pi*C_SI.c/lambda;
% Up2=9.337287*peak_intensity*1e-18*(lambda*1e6)^2
Up = C_SI.qE^2/(4*C_SI.mE)*(FE*field)^2/omega^2/C_SI.qE; % [eV]
E_cutoff = 10.007*Up + 0.538*WorkFunction; % [eV]
% Keldysh parameter
% VI=10*qE; % typical 10 eV ionization energy in Joules
% VI = 5*C_SI.qE; % 5 eV work function in Joules
keldysh = omega*sqrt(2*C_SI.mE*WorkFunction*C_SI.qE)/(C_SI.qE*FE*field);
clc;
disp(['Peak intensity = ' num2str(peak_intensity*1e-4,'%1.3e') ' W/cm^2']);
disp(['Peak field = ' num2str(field,'%1.3e') ' V/m']);
disp(['Up = ' num2str(Up,'%1.3f') ' eV']);
disp(['Keldysh gamma = ' num2str(keldysh,'%1.3f')]);
disp(['E_{cutoff} = ' num2str(E_cutoff,'%1.3f') ' eV']);
% This code simulates the Phoibos 100 Hemispherical Energy Analyzer (HEA).
% The analysis is based on [https://en.wikipedia.org/wiki/Kepler_orbit#math_1]
% see section 'Determination of the Kepler orbit that corresponds to a given initial state'
constants_si;
k1 = 1/4/pi/C_SI.eps0;
R0 = 0.1; % [m] 100 mm for Phoibos 100
R_inner = 0.75*R0;
R_outer = 1.25*R0;
E_kin_eV = 1e-8; % [eV]
E_kin = E_kin_eV * C_SI.qE;
v_mag = sqrt(E_kin*2/C_SI.mE);
gamma = pi/2; % angle relative to radial direction in [rad]
E_pass = 50; % pass energy in [eV]
Delta_U = E_pass; % potential difference between inner and outer spheres.
q0 = Delta_U/k1/(1/R_inner - 1/R_outer); % Charge on inner plate
kappa = q0*k1;

r_0 = R0;
v_r = v_mag * cos(gamma);
v_t = v_mag * sin(gamma);
p = (r_0*v_t)^2/kappa;
v_0 = sqrt(kappa/p);
phi_0 = atan(v_r/(v_t-v_0));
excentricity = v_r/v_0/sin(phi_0);
% r_0b = p./(1+excentricity*cos(phi0));

if excentricity == 0
    disp('Orbit type: circular');
elseif excentricity>0 && excentricity<1
    disp('Orbit type: elliptical');
elseif excentricity == 1
    disp('Orbit type: parabolic');
elseif excentricity >= 1
    disp('Orbit type: hyperbolic');
end

phi_vec = 0:2*pi/99:2*pi;
r_vec = p./(1+excentricity*cos(phi_vec));
r_det1 = p./(1+excentricity*cos(phi_0+pi));
x_vec = r_vec.*cos(phi_vec);
y_vec = r_vec.*sin(phi_vec);
x_0 = r_0.*cos(phi_0);
y_0 = r_0.*sin(phi_0);
% figure;plot(phi_vec,r_vec,'ko');hold on;plot(phi_0,r_0,'ro')
figure;plot(x_vec,y_vec,'ko');hold on;plot(x_0,y_0,'ro')
%%
r_det2 = (E_kin - E_pass)/E_pass*2*R0 + R0;
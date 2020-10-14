function out = calc_dscan(varargin)

omega = varargin{1}; % [rad/fs]
E_om = varargin{2}; % Electric field vs omega [arb. units]
switch nargin
    case 2
        material = 'FS';
        mat_vec = [-3:0.1:3];
    case 3
        material = varargin{3};
        mat_vec = [-3:0.1:3];
    case 4
        material = varargin{3};
        mat_vec = varargin{4};
end
N_mat = length(mat_vec);

wl = map2colvec(2*pi*300./omega); % [nm]
n_mat = sellmeier(wl*1e-3,material);
Phi = 2*pi./(wl*1e-6).*n_mat; % spectral phase vs wavelength for unit thickness of material
Phi = Phi - (Phi(end)-Phi(1))/(omega(end)-omega(1))*omega;
Phi = Phi - min(Phi);

crystal_length = 0*12; % [um]
crystal_theta = 27.5; % [deg]

for ind1 = 1:N_mat
    Phi_total = mat_vec(ind1)*Phi;
    [E_omega2,I_lambda] = calc_SHG_spectrum([map2colvec(omega) map2colvec(abs(E_om)) map2colvec(unwrap(angle(E_om.*exp(1i*Phi_total))))],crystal_length,crystal_theta);
    trace_Dscan_lambda(ind1,:) = map2rowvec(I_lambda(:,2));
    trace_Dscan_omega(ind1,:) = map2rowvec(E_omega2(:,2));
%     disp([num2str(ind1) '/' num2str(N_mat)]);
end

out.trace_lambda = trace_Dscan_lambda;
out.trace_omega = trace_Dscan_omega;
out.axis_lambda = map2rowvec(I_lambda(:,1));
out.axis_omega = map2rowvec(E_omega2(:,1));
out.mat_vec = mat_vec;
end
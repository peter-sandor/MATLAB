function main_NLI(nu,n,amp0)
% Solving coupled wave equations for SFG
% Source: W. R. Boyd Nonlinear Optics (3rd ed. 2007)

% Input arguments:  nu - [3x1] vector, frequencies for the fields in [PHz]
%                   amp0 - [3x1] vector, relative initial amplitudes for the fields
%                   [V/nm]
%                   n - [3x1] refractive indices

c=300; % [nm/fs]
z_max = 5e3; % [nm]
z_min = 0;
N = 1*256;
% z_stp = (z_max-z_min)/(N-1);
% z = (z_min:z_stp:z_max)';
k=(2*pi)*map2colvec(nu).*map2colvec(n)/c; % wave vector magnitudes (k=n*omega/c)
delta_k=k(1)+k(2)-k(3); % Boyd 2.2.11
d_eff=2e-2; % [nm/V] 
param=2*d_eff/c^2*((2*pi)*map2colvec(nu)).^2./k;
% amp0 = [1;0.2;0];

% from Boyd Nonlinear Optics book, 3rd edition, equations 2.2.10, 2.2.12a
% and 2.2.12b
func_NLI=@(z,A) [1i*param(1)*A(3)*conj(A(2))*exp(-1i*delta_k*z); 1i*param(2)*A(3)*conj(A(1))*exp(-1i*delta_k*z); 1i*param(3)*A(1)*A(2)*exp(+1i*delta_k*z);];

[z, amplitudes] = ode45(func_NLI,[z_min z_max],amp0);
intensity =1/2*8.85e-8*3*(1e9*abs(sum(amplitudes,2))).^2; % total intensity in [W/cm^2]

% subplot(211)
plot(z, abs(amplitudes(:,1)).^2, 'r-')
hold on
plot(z, abs(amplitudes(:,2)).^2, 'r--')
plot(z, abs(amplitudes(:,3)).^2, 'b-')
hold off
xlabel('x [nm]')
ylabel('field amplitude [V/nm]')
legend(num2str(map2colvec(nu)))
title('field amplitudes vs length; legend: frequencies in [PHz]')
% subplot(212)
% plot(z,intensity,'k-')
% xlabel('x [nm]')
% ylabel('Total Intensity [W/cm^2]')
end
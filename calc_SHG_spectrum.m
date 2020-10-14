function [varargout] = calc_SHG_spectrum(E_om,crystal_length,crystal_theta)
% This code calculates SHG/SFG spectrum (field magnitude and phase vs omega OR intensity vs wavelength)
% for type I phasematching if fundamental excitation (field magnitude and phase vs. omega) is known.
% Pump depletion is not included.
% Input variables:      - E_om is an Nx3 vector, with three columns:
%                           first: abscissa vector omega
%                           second: electric field amplitude vs omega
%                           third: electric field phase vs omega
%                       - crystal_length in microns
%                       - crystal_angle in degrees
% Output variables (optional):
%                       - Nx3 array containing the following columns:
%                           first: abscissa vector omega [rad/fs]
%                           second: electric field amplitude vs omega
%                               (taking into account finite phasematching bandwidth)
%                           third: electric field amplitude vs omega
%                               (assuming infinite phasematching bandwidth)
%                       - Nx3 array containing the following columns:
%                           first: abscissa vector lambda [nm]
%                           second: intensity vs lambda
%                               (taking into account finite phasematching bandwidth)
%                           third: intensity vs lambda
%                               (assuming infinite phasematching bandwidth)

%% 
[t1,E1,omega,E_omega] = apply_spectral_mask_optimizing(E_om(:,1),E_om(:,2).*exp(1i*E_om(:,3)),'1');

if crystal_length~=0 && crystal_theta~=0 % Calculate using autoconvolution of fundamental spectrum
    source = 2;
    crystal_theta = deg2rad(crystal_theta); % Convert angle to degrees
    N0 = length(omega);
    om_step = omega(2)-omega(1);
    max_index = round(max(vec2ind(abs(E_omega)>0))*1.2);
    % ind_select0 = (N0/2-max_index:N0/2+max_index);
    ind_select0 = [1:max_index N0-max_index+1:N0];
    omega1 = fftshift(omega(ind_select0));
    E_omega1 = fftshift(E_omega(ind_select0));
    N = length(omega1);
    lambda1 = 300./omega1*2*pi*1e-9;
    omega1_flip = flip(omega1,1);
    omega1_flip = [max(omega1_flip)+om_step; omega1_flip(1:end-1)];
    E_omega1_flip = flip(E_omega1,1);
    E_omega1_flip = [0; E_omega1_flip(1:end-1)];
    ind_check = vec2ind(abs(E_omega1) == max(abs(E_omega1)));

    for ind1 = 1:N
    %     omega2(ind1) = (ind1 - 1)*om_step;
        omega1_shift = circshift(omega1_flip,[ind1-1 0]);
        E_shift = circshift(E_omega1_flip,[ind1-1 0]);
        omega_SHG(ind1) = omega1(ind_check) + omega1_shift(ind_check);
        lambda_shift = 300./omega1_shift*2*pi*1e-9;
        lambda_SHG(ind1) = 300./omega_SHG(ind1)*2*pi*1e-9;
        % Refractive index of extraordinary axis as a function of the angle theta
        n_e_theta(ind1) = (n_e(lambda_SHG(ind1),source).*n_o(lambda_SHG(ind1),source))./(sqrt(n_o(lambda_SHG(ind1),source).^2.*sin(crystal_theta).^2 + n_e(lambda_SHG(ind1),source).^2.*cos(crystal_theta).^2));
        % Phase-mismatch
        Delta_k(:,ind1) = map2colvec((2*pi./lambda1.*n_o(lambda1,source) + (2*pi./lambda_shift.* n_o(lambda_shift,source)) - 2*pi./(lambda_SHG(ind1))*n_e_theta(ind1)));
        p_m(:,ind1) = sin(1/2*Delta_k(:,ind1).*crystal_length*1e-6).^2./(1/2*Delta_k(:,ind1)*crystal_length*1e-6).^2;
        ind_select(:,ind1) = ~any([isnan(p_m(:,ind1)) isinf(p_m(:,ind1))],2);
        if max(double(ind_select(:,ind1)))
            % calculate case for imperfect phase matching due to crystal
            % bandwidth
            E_SHG(ind1) = sum(p_m(ind_select(:,ind1),ind1).*E_omega1(ind_select(:,ind1)).*E_shift(ind_select(:,ind1)));
        else
            E_SHG(ind1) = 0;
        end
    end
else % Calculate using Fourier transform of the time-dependence of field squared (infinite phase-matching bandwidth)
    E_SHG = map2colvec(fft(E1.^2));
    I_SHG = abs(E_SHG).^2;
    omega_SHG = omega;
end

% ind_select2 = vec2ind(E_SHG~=0);
% omega_SHG = omega_SHG(ind_select2);
% E_SHG = E_SHG(ind_select2);
% omega_range = 2*[2 3];
% index2 = all([map2colvec(omega_SHG>omega_range(1)) map2colvec(omega_SHG<omega_range(2))],2);
index2 = vec2ind(omega_SHG>0);
omega_plot = omega_SHG(index2);
lambda_plot = 2*pi*300./(omega_plot);
% norm_factor0 = max(abs(E_SHG0(index2)).^2.*omega_plot.^2);
norm_factor0 = 1;
% I_SHG_plot_0 = abs(E_SHG0(index2)).^2.*omega_plot_2.^2/norm_factor0;
I_SHG_plot = abs(E_SHG(index2)).^2.*omega_plot.^2/norm_factor0;

if 0
    figure;plot(lambda_plot,I_SHG_plot,'k')
end

varargout{1} = [map2colvec(omega_SHG) map2colvec(E_SHG)];
varargout{2} = [map2colvec(lambda_plot) map2colvec(I_SHG_plot)];
end

% Source 2: Eimerl et al. 1987: BBO; n(o) 0.22-1.06 µm, https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Eimerl-o
% Source 3: Zhang et al. 2000: BBO; n(o) 0.64-3.18 µm, https://refractiveindex.info/?shelf=main&book=BaB2O4&page=Zhang-o
function n = n_o(lambda,source)
% N_BBO_O Refractive index of o-wave in BBO.
lambda = lambda * 1e6;
if source == 1
   n = sqrt(2.7359 + 0.01878 ./ (lambda.^2 - 0.01822) - 0.01354 * lambda.^2);
elseif source == 2
   n = sqrt(2.7405 + 0.0184 ./ (lambda.^2 - 0.0179) - 0.0155 * lambda.^2);
   n(any([map2colvec(lambda<0.22) map2colvec(lambda>1.06)],2)) = NaN;
elseif source == 3
   n = sqrt(2.7359 + 0.01878 ./ (lambda.^2 - 0.01822) ...
      - 0.01471 * lambda.^2 + 0.0006081*lambda.^4 - 0.00006740*lambda.^6);
end
end

function n = n_e(lambda,source)
% N_BBO_O Refractive index of e-wave in BBO. Lambda in m.
lambda = lambda * 1e6;
lambda(lambda<0) = NaN;
if source == 1
   n = sqrt(2.3753 + 0.01224 ./ (lambda.^2 - 0.01667) - 0.01516 * lambda.^2);
elseif source == 2
   n = sqrt(2.3730 + 0.0128 ./ (lambda.^2 - 0.0156) - 0.0044 * lambda.^2);
   n(any([map2colvec(lambda<0.22) map2colvec(lambda>1.06)],2)) = NaN;
elseif source == 3
   n = sqrt(2.3753 + 0.01224 ./ (lambda.^2 - 0.01667) ...
      - 0.01627 * lambda.^2 + 0.0005716*lambda.^4 - 0.00006305*lambda.^6);
end
end

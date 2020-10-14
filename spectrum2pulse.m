function varargout = spectrum2pulse(omega_in,field_in,M)
% omega_in:   angular frequencies
% field_in:   field in spectral domain (can be complex)
% M:        string describing the field mask in the spectral domain as a function of frequency. Use
%           'omega' to refer to angular frequency as the independent variable and use omega0 to
%           denote the central angular frequency (to be determined in this code).

N_iter = 4;
SampleOm = 64;
N0=size(field_in,1);
omega=omega_in;
omega_min0=min(omega_in);
omega_max0=max(omega_in);
omega_range0=omega_max0-omega_min0;
omega_step0=omega_range0/(N0-1);
omega0 = sum(omega_in.*abs(field_in).^2)/sum(abs(field_in).^2);

[mean_om, sigma_om, FWHM_om, ind_FWHM, ind_crop] = peak_props([omega_in (abs(field_in)/max(abs(field_in))).^2]);

T_period = 2*pi/mean_om; % [fs]
points_per_cycle = 10;
delta_t = T_period/points_per_cycle;
tmax = 2*pi*SampleOm/FWHM_om;
Nt = 2^(ceil(log2(floor(2*tmax/delta_t))));
ind_iter = 1;
ind_nonzero = [1 Nt];
while (ind_nonzero(2)-ind_nonzero(1)>round(Nt*3/4)) && (ind_iter<=N_iter)
    if ind_iter > 1
        tmax = tmax*2;
    end
    Nt = 2^(ceil(log2(floor(2*tmax/delta_t))));
    taxis = map2colvec(-tmax:2*tmax/(Nt-1):tmax);
    omega = map2colvec(FourierAxis(taxis));
    sp_omega = zeros(size(omega));
    index_insert1 = min(vec2ind(omega(1:Nt/2)>omega_min0));
    index_insert2 = max(vec2ind(omega(1:Nt/2)<omega_max0));
    % Use 'spline' interpolation for complex numbers, or interpolate real and imaginary parts separately.
    sp_omega(index_insert1:index_insert2) = interp1(omega_in,field_in,omega(index_insert1:index_insert2),'spline');
    spfield = apply_spectral_mask(omega,sp_omega,M);
%     eval(['mask=' M ';']);
%     spfield0 = sp_omega.*map2colvec(mask);
%     spfield = spfield0;
    % Fourier-transform to time domain
    tdfield = ifftshift(ifft(spfield));
    [~, ~, ~, ~, ind_nonzero] = peak_props([taxis (abs(tdfield)/max(abs(tdfield))).^1]);
    ind_iter = ind_iter+1;
    if ind_iter == N_iter+1
        disp('Check time delay range of time-domain waveform.')
    end
end

if nargout>0
    varargout{1} = map2colvec(taxis);
    varargout{2} = map2colvec(tdfield);
    varargout{3} = map2colvec(omega);
    varargout{4} = map2colvec(spfield);
    varargout{5} = omega0;
end
end


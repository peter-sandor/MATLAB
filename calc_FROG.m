function [spectrax,tauax,signal_FROG] = calc_FROG(timeax,field0,type,freq_or_lambda,varargin)

% This code calculates SHG, SD and interferometric SHG FROG traces from the
% time dependent electric field 'field0'. The resulting trace can be either
% as a function of frequency (in THz) or in wavelength (in nm).
% Input parameters:     - timeax: abscissa vector (in fs) for the electric field (Nx1 real array)
%                       - field0: Nx1 complex electric field as a function
%                       of time (Nx1 complex array)
%                       - type: string specifying the signal to calculate
%                       ('SHG', 'SD' or 'IFROG')
%                       - freq_or_lambda: string containing either 'omega'
%                       or 'lambda' to switch between spectral axes
%                       - time axis for interpolation (Mx1 array)
%                       - spectral axis for interpolation (Mx1 array)
% Output parameters:    - spectrax: spectral axis (either angular frequency in rad/fs or wavelength in nm)
%                       - tauax: time axis in 'fs'
%                       - signal_SHG: noncollinear SHG signal (MxM real array)
%                       - signal_SD: noncollinear SD signal (MxM real array)
%                       - signal_IFROG: collinear (interferometric) FROG signal (MxM real array)

c = 300; % [nm/fs]
Nt = 2^ceil(log2(length(timeax)));
timeax = timeax - mean(timeax(abs(field0)==max(abs(field0))));
timeax_range = timeax(end)-timeax(1);
tauax_step = timeax_range/(Nt-1);
tauax = map2colvec(-Nt/2*tauax_step:tauax_step:(Nt/2-1)*tauax_step);
field1 = map2colvec(interp1(timeax,field0,tauax,'linear'));
field1(isnan(field1)) = 0;
omega = map2colvec(FourierAxis(tauax));
% field1_fft = fft(field1);
% omega_center = peak_props([map2colvec(omega(1:Nt/2)) map2colvec(abs(field1_fft(1:Nt/2)).^2)]);

field2 = [zeros([0*Nt/2 1]); field1; zeros([0*Nt/2 1])];
signal_FROG = zeros([Nt Nt]);

for ind1=1:Nt
    field2_shifted = circshift(field2,[Nt/2-ind1 1]);
    switch type
        case 'IFROG'
            signal_FROG(:,ind1) = map2rowvec(abs(fft((field1 + field2_shifted).^2)).^2);
        case 'SHG'
            signal_FROG(:,ind1) = map2rowvec(abs(fft(field1.*field2_shifted)).^2);
        case 'SD'
            signal_FROG(:,ind1) = map2rowvec(abs(fft(field1.^2.*conj(field2_shifted))).^2);
        otherwise
            signal_FROG(:,ind1) = [];
    end
end
if strcmp(freq_or_lambda,'omega')==1
    spectrax = omega;
    signal_FROG=signal_FROG/max(max(signal_FROG));
elseif strcmp(freq_or_lambda,'lambda')==1
    spectrax = flip(2*pi*c./omega(1:round(Nt/2)),1); % spectrax = lambda
    spectrax(spectrax == Inf) = NaN;
    omega_sq_factor = repmat(omega.^2,[1 Nt]);
    signal_FROG = signal_FROG.*omega_sq_factor; % apply Jacobian: f(lambda) = g(omega)*omega^2;
    signal_FROG = flip(signal_FROG(1:round(Nt/2),:)/max(max(signal_FROG(1:round(Nt/2),:))),1);
end
if nargin == 6
    tauax_ip = varargin{1};
    spectrax_ip = varargin{2};
    [TAU_orig, SP_orig] = meshgrid(tauax,spectrax(1:round(Nt/2)));
    [TAU_IP, SP_IP] = meshgrid(tauax_ip,spectrax_ip);
    tauax = tauax_ip;
    spectrax = spectrax_ip;
    signal_FROG = interp2(TAU_orig,SP_orig,signal_FROG(1:round(Nt/2),:),TAU_IP,SP_IP,'linear',0);
end
end
function varargout = simulate_pulse_shaper(datain,M)

% datain:   first column --> omega
%           second column --> field in spectral domain
% M:        string describing the field mask in the spectral domain as a function of frequency. Use
%           'omega' to refer to different angular frequencies and use omega0 to
%           denote the central angular frequency.

SampleT=10;
SampleOm=40;
N0=size(datain,1);
omega=datain(:,1);
omega_min0=min(datain(:,1));
omega_max0=max(datain(:,1));
omega_range0=omega_max0-omega_min0;
omega_step0=omega_range0/(N0-1);
omega0=sum(datain(:,1).*datain(:,2))/sum(datain(:,2));

[mean_om sigma_om FWHM_om ind_FWHM ind_crop] = calc_stats([datain(:,1) (abs(datain(:,2))/max(abs(datain(:,2)))).^2]);

delta_t=2*pi/(omega0*SampleT);
tmax=2/(FWHM_om/SampleOm);
N1=floor(3*tmax/delta_t);
Nt=2^(ceil(log2(N1)));
t=map2colvec(-tmax:2*tmax/(Nt-1):tmax);
omega_step=omega0*SampleT/N1;
omega_max=omega_step*(Nt/2-1);
omega_min=-omega_step*(Nt/2);
omega=map2colvec(0:omega_step:omega_max);
sp_omega=zeros([Nt/2 1]);
index_insert1=min(vec2ind(omega>omega_min0));
index_insert2=max(vec2ind(omega<omega_max0));
sp_omega(index_insert1:index_insert2)=interp1(datain(:,1),datain(:,2),omega(index_insert1:index_insert2),'linear');

eval(['mask=' M ';']);
spfield0=sp_omega.*map2colvec(mask);
% spfield=[0; flipdim(spfield0(2:end),1); spfield0];
spfield=[zeros(size(spfield0)); spfield0];
omega=[map2colvec(omega_min:omega_step:-omega_step); omega];

taxis=FourierAxis(omega);
% taxis=taxis-mean(taxis);
% Fourier-transform to time domain
tdfield=ifftshift(ifft(fftshift(spfield)));

if nargout>0
    varargout{1}=map2colvec(ifftshift(taxis));
    varargout{2}=map2colvec(tdfield);
    varargout{3}=map2colvec(omega);
    varargout{4}=map2colvec(spfield);
end
end


function spectrum2pulse(varargin)

% code that calculates the pulse intensity in the time domain given an
% experimentally measured spectra (spectral density vs. wavelength) and
% assumed phase coefficients of different orders in the frequency domain
% usage:
% spectrum2pulse -  calls spproc to open spectrum file and load data;
%                   assumes flat phase throughout the spectrum
% spectrum2pulse(data) -    the argument 'data' contains the wavelengths in
%                           the first column and intensities in the second
% spectrum2pulse(data,phasecoeffs) -    phasecoeffs is a vector, where the
%                                       n-th element is the n-th order phase coefficient d^n(Phi)/d(omega)^n

if nargin==0
    data=spproc;
elseif nargin>=1
    data=varargin{1};
end
[data2(:,1),data2(:,2)]=consolidator(data(:,1),data(:,2),'mean');
omega1=2*pi*300./data2(:,1); % if wavelength is in [nm], omega is in [rad/fs]
omega_min=min(omega1);
omega_max=max(omega1);
omega_span=omega_max-omega_min;
N=size(data2,1);
% spint=data(:,2);
if N>128
    N2=128;
elseif N<=128
    N2=N;
end
omega_step=omega_span/(N2-1);
omega=map2colvec(0:omega_step:omega_span);
spint=interp1(omega1-omega_min,data2(:,2),omega,'linear');
spfield=sqrt(abs(spint));
zeropadfactor=4;
if nargin==2
    omega_mean=sum(spfield.*omega)/sum(spfield);
    spfield=spfield.*exp(+i*pi*double(omega>=omega_mean+varargin{2}*omega_step));
end
% figure;plot(omega,abs(spfield),'k-')
% hold on;plot(omega,unwrap(angle(spfield))/2/pi,'r-')
% hold off
spfield=[spfield; zeros([(zeropadfactor-1)*N2 1])];

tdfield=ifftshift(ifft(spfield));
taxis=2*pi*map2colvec(0:1/(zeropadfactor*omega_span):(zeropadfactor*N2-1)/(zeropadfactor*omega_span));
taxis=taxis-mean(taxis);
hndl2=figure;plot((abs(tdfield)/max(abs(tdfield))).^2,'k','linewidth',1)
title('Select region of interest')
temp=ginput(2);
crop_start=uint16(round(min(temp(1,1),temp(2,1))));
crop_end=uint16(round(max(temp(1,1),temp(2,1))));
plot(taxis(crop_start:crop_end),(abs(tdfield(crop_start:crop_end))/max(abs(tdfield(crop_start:crop_end)))).^2,'k','linewidth',2)
hold on;plot(taxis(crop_start:crop_end),0.5*ones(size(taxis(crop_start:crop_end))),'r')
xlabel('Delay [fs]')
title('Choose points for FWHM')
temp=ginput(2);
title('')
% close(hndl1);
FWHM=abs(diff(temp(:,1)));
disp(['temporal FWHM = ' num2str(FWHM) ' fs'])



end
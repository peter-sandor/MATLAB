function varargout = spectrum2pulse_old(varargin)

% code that calculates the pulse intensity in the time domain given an
% experimentally measured spectra (spectral density vs. wavelength) and
% assumed phase coefficients of different orders in the frequency domain
% usage:
% spectrum2pulse -  calls spproc to open spectrum file and load data;
%                   assumes flat phase throughout the spectrum
% spectrum2pulse(data) -    the argument 'data' contains the wavelengths in
%                           the first column and intensities in the second
%                           if phases are also available, that can be input two ways:
%                            i) phase goes into 3rd column
%                           ii) the real part of the spectrum goes in the 3rd, imaginary in the 4th column
% spectrum2pulse(data,phasecoeffs) -    phasecoeffs is a vector, where the
%                                       n-th element is the n-th order phase coefficient d^n(Phi)/d(omega)^n
%                                       this phase will be applied in addition to what was supplied (if any) in
%                                       'data'.

if nargin==0
    data=spproc;
elseif nargin>=1
    data=varargin{1};
end
[data2(:,1),data2(:,2:max(2,size(data,2)))]=consolidator(data(:,1),data(:,2:max(2,size(data,2))),'mean');
data2(:,2)=data2(:,2).*data2(:,1).^2; % apply Jacobian: g(omega)= lambda^2*f(lambda)
data2(:,2)=data2(:,2)/max(data2(:,2));
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
omega=map2colvec(omega_min:omega_step:omega_max);
data_ip=interp1(omega1,data2(:,2:max(2,size(data,2))),omega,'linear');
data_ip(:,1)=sqrt(abs(data_ip(:,1))); % convert from intensity to field strength
if size(data_ip,2)==2
    spfield=data_ip(:,1).*exp(i*data_ip(:,2));
elseif size(data_ip,2)==3
    spfield=data_ip(:,1).*exp(i*(unwrap(angle(data_ip(:,2)+i*data_ip(:,3)))));
else 
    spfield=data_ip(:,1);
end
zeropadfactor=4;
if nargin==2
    omega_mean=sum(abs(spfield).*omega)/sum(abs(spfield));
    for ind1=1:length(varargin{2})
        spfield=spfield.*exp(-i/factorial(ind1)*varargin{2}(ind1)*(omega-omega_mean).^ind1);
    end
end
% figure;plot(omega,abs(spfield),'k-')
% hold on;plot(omega,unwrap(angle(spfield))/2/pi,'r-')
% hold off
spfield=[spfield; zeros([(zeropadfactor-1)*N2 size(spfield,2)])];

tdfield=ifftshift(ifft(spfield));
taxis=2*pi*map2colvec(0:1/(zeropadfactor*omega_span):(zeropadfactor*N2-1)/(zeropadfactor*omega_span));
taxis=taxis-mean(taxis);
hndl2=figure;plot((abs(tdfield)/max(abs(tdfield))).^2,'k','linewidth',1)
title('Select region of interest')
temp=ginput(2);
crop_start=uint16(round(min(temp(1,1),temp(2,1))));
crop_end=uint16(round(max(temp(1,1),temp(2,1))));
close(hndl2)

figure;
plot(taxis(crop_start:crop_end),(abs(tdfield(crop_start:crop_end))/max(abs(tdfield(crop_start:crop_end)))).^2,'k','linewidth',2)
hold on;plot(taxis(crop_start:crop_end),0.5*ones(size(taxis(crop_start:crop_end))),'r')
xlabel('Delay [fs]')
title('Choose points for FWHM')
temp=ginput(2);
title('')
% close(hndl1);
FWHM=abs(diff(temp(:,1)));
disp(['temporal FWHM = ' num2str(FWHM) ' fs'])

if nargout>0
    taxis_out=taxis(crop_start:crop_end);
    fieldout=tdfield(crop_start:crop_end)/max(tdfield(crop_start:crop_end));
    Nout=length(taxis_out);
    taxis_ip=taxis_out(1):(taxis_out(end)-taxis_out(1))/127:taxis_out(end);
    field_ip=interp1(taxis_out,fieldout,taxis_ip);
    varargout{1}=[map2colvec(taxis_ip) map2colvec(field_ip)];
end

end
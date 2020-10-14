function varargout = SpectrumConvert(varargin)

% code that converts given the intensity measured as a function of wavelength (in nanometers) to intensity as a function of natural frequency (omega)
% it also converts (resamples) the spectral phase if it is passed along the
% spectral densities

% SpectrumConvert -  calls spproc2 to open spectrum file and load data;
%                   assumes flat phase throughout the spectrum
% SpectrumConvert(data) -    the argument 'data' contains the wavelengths in
%                           the first column and intensities in the second
%                           if phases are also available, that can be included as a 3rd column
% SpectrumConvert(data,N) -    the number of points 'N', used for
%                           interpolation, is also specified (default: N=128)
% out = SpectrumConvert - the variable out contains the natural frequencies
%                           in the first column, electric field strengths
%                           in the second, and (if available) phases in the
%                           third
if nargin==0
    data=spproc2;
    N2=2*256;
    flag_norm=0;
elseif nargin==1
    data=varargin{1};
    N2=2*256;
    flag_norm=0;
elseif nargin==2
    data=varargin{1};
    N2=varargin{2};
    flag_norm=0;
elseif nargin==3
    data=varargin{1};
    N2=varargin{2};
    flag_norm=varargin{3};
end
thrs=0.02;
[data2(:,1),data2(:,2:max(2,size(data,2)))]=consolidator(data(:,1),data(:,2:max(2,size(data,2))),'mean');
% data2(:,2)=data(:,2)-mean(data(1:10,2),1);
data2(:,2)=data2(:,2).*data2(:,1).^2; % apply Jacobian: g(omega)= lambda^2*f(lambda)
if flag_norm
    data2(:,2)=data2(:,2)/max(data2(:,2));
end
omega1=2*pi*300./data2(:,1); % if wavelength is in [nm], omega is in [rad/fs]
omega_min=min(omega1);
omega_max=max(omega1);
omega_span=omega_max-omega_min;
N=size(data2,1);
% if N<128
%     N2=128;
% elseif N>=128
%     maxpower=10;
%     temp1=abs(2.^([0:maxpower])-N);
%     N2=2^max(unique((temp1==min(temp1)).*(0:maxpower)));
% end
omega_step=omega_span/(N2-1);
omega_ip=map2colvec(omega_min:omega_step:omega_max);
data_ip=interp1(omega1,data2(:,2:max(2,size(data,2))),omega_ip,'linear');
data_ip(data_ip(:,1)<thrs,1)=0;
data_ip(:,1)=sqrt(abs(data_ip(:,1))); % convert from intensity to field strength
varargout{1}=[map2colvec(omega_ip) data_ip(:,1:max(1,size(data_ip,2)))]; % omega and E(omega)
varargout{2}=[data(:,1) sqrt(abs(data(:,2)))]; % lambda and E(lambda) (not intensity, but field!)
end
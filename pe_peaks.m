function varargout = pe_peaks(varargin)

% [filename,pathname,filterindex] = uigetfile;
% data=load(strcat(pathname,filename));
if nargin==0
    [FileName,PathName] = uigetfile({'*.sco*','OceanOptics spectrum file';'*.txt','Text File';'*.*','All Files (*.*)'},'Pick a spectrum file');
    S=importdata(strcat(PathName,FileName));
    % S=uiimport('-file',{'*.scope','OceanOptics spectrum file';'*.txt','text data';'*.*','All Files (*.*)'},'Pick a spectrum file');
    if isstruct(S)==1
        data=S.data;
    elseif isnumeric(S)==1
        data=S;
    end
    clear S;
    hndl1=figure;plot(data(:,1),data(:,2),'k')
    title('Choose start and end of background')
    temp=ginput(2);
    bckg_start=min(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
    bckg_end=max(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
    if bckg_start < 1
        bckg_start=1;
    end
    if bckg_end > size(data,1)
        bckg_end=size(data,1);
    end        
%     bckg_start=min(uint16(round(temp(:,1))));
%     bckg_end=max(uint16(round(temp(:,1))));
    close(hndl1)
    photon_nrg=1240000/meanwl_Ver2(data,bckg_start,bckg_end) % [meV]
else photon_nrg = varargin{1};
end
temp=ginput(2);
calib_parameter=abs(photon_nrg/(temp(1,1)^2-temp(2,1)^2)); % [meV]
varargout{1}=calib_parameter;
varargout{2}=(temp(1,1)+temp(2,1))/2;
end
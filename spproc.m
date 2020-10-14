function varargout = spproc(varargin)
% Code for reading in spectrum files, cropping the data and determining
% FWHM by eyeballing (ginput)
% Can handle Oceanoptics spectrum files (.sco and Master.Scope) and
% ASCII-delimited text (.txt)
% Filename and full path should be supplied in argument, if not, dialog
% window pops up. Output (optional) is an Mx2 double array (first column: wavelengths,
% second column: corresponding spectral densities)
if nargin==1
%     [PathName, name, ext, versn] = fileparts(varargin{1});
%     PathName=strcat(PathName,'\');
%     FileName=strcat(name,ext);
    if ischar(varargin{1})
        S=importdata(varargin{1});
    elseif ishandle(varargin{1})
        hndl_in=varargin{1};
        [FileName,PathName] = uigetfile({'*.*','All Files (*.*)';'*.txt','Text File';},'Pick a spectrum file');
        S=importdata(strcat(PathName,FileName));
    else 'Wrong type of argument.'
        return;
    end
elseif nargin==0
    [FileName,PathName] = uigetfile({'*.*','All Files (*.*)';'*.txt','Text File';},'Pick a spectrum file');
    S=importdata(strcat(PathName,FileName));
end
% S=importdata(strcat(PathName,FileName));
% S=uiimport('-file',{'*.scope','OceanOptics spectrum file';'*.txt','text data';'*.*','All Files (*.*)'},'Pick a spectrum file');
if isstruct(S)==1
    temp=S.data;
elseif isnumeric(S)==1
    temp=S;
end
[data(:,1),data(:,2)]=consolidator(temp(:,1),temp(:,2),'mean');
clear S temp;
hndl1=figure;plot(data(:,2),'k')
title('Choose start and end of background')
temp=ginput(2);
bckg_start=uint16(round(min(temp(1,1),temp(2,1))));
bckg_end=uint16(round(max(temp(1,1),temp(2,1))));
% bckg_start=min(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
% bckg_end=max(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
if bckg_start < 1
    bckg_start=1;
end
if bckg_end > size(data,1)
    bckg_end=size(data,1);
end        
title('Select region of interest')
temp=ginput(2);
crop_start=uint16(round(min(temp(1,1),temp(2,1))));
crop_end=uint16(round(max(temp(1,1),temp(2,1))));
if crop_start < 1
    crop_start=1;
end
if crop_end > size(data,1)
    crop_end=size(data,1);
end      
spectrum=(data(crop_start:crop_end,2)-mean(data(bckg_start:bckg_end,2)))/max(data(crop_start:crop_end,2)-mean(data(bckg_start:bckg_end,2)));
wl=data(crop_start:crop_end,1);
c=300; % [nm/fs]
% meanwl=sum(wl.*spectrum)/sum(spectrum);
% std_dev=std(hist(spectrum,wl-meanwl));
% plot(wl,spectrum,'k.-')
% hold on;plot(wl,0.5*ones(size(spectrum)),'r')
% hold off;
% title('Choose points for FWHM')
% temp=ginput(2);
% title('')
close(hndl1);
% FWHM=abs(diff(temp(:,1)));
[wlmean std_dev FWHM]=calc_stats([map2colvec(wl) map2colvec(spectrum)]);
meanfreq=c/wlmean;
FWHM2=300/wlmean^2*FWHM;
disp(['mean wavelength = ' num2str(wlmean) ' nm (= ' num2str(meanfreq) ' PHz)'])
disp(['FWHM = ' num2str(FWHM) ' nm (= ' num2str(FWHM2) ' PHz)'])
if exist('hndl_in')==1
    answer = inputdlg('Enter a line color identifier');
    hold on;plot(hndl_in,wl,spectrum,[answer{1} '.-'])
    hold off
end
if nargout==1
    varargout{1}=[map2colvec(wl),map2colvec(spectrum)];
%     close(hndl1);
end
end
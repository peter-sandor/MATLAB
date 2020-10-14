function varargout = read_spectrum(varargin)
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
        file_path = varargin{1};
    else 'Wrong type of argument.'
        return;
    end
elseif nargin==0
    [FileName,PathName] = uigetfile({'*.*','All Files (*.*)';'*.txt','Text File';},'Pick a spectrum file');
    file_path = [PathName FileName];
end
S=importdata(file_path);
% S=uiimport('-file',{'*.scope','OceanOptics spectrum file';'*.txt','text data';'*.*','All Files (*.*)'},'Pick a spectrum file');
if isstruct(S)==1
    temp=S.data;
elseif isnumeric(S)==1
    temp=S;
end
[data(:,1),data(:,2)]=consolidator(temp(:,1),temp(:,2),'mean');

if nargout == 1
    varargout{1}=data;
elseif nargout == 0
    figure;plot(data(:,1),data(:,2),'k')
    xlim([min(data(:,1)) max(data(:,1))]);
    xlabel('wavelength [nm]')
    setfigP;
end
end
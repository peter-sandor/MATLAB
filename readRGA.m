function [dataout,varargout] = readRGA(varargin)

if nargin==0
    [filename,pathname]=uigetfile({'*.xml','XML files';'*.*','All files'});
    full_path=strcat(pathname,filename);
elseif nargin==1
    full_path=varargin{1};
end
fileid=fopen(full_path);
header=textscan(fileid,'<Data LowMass="%f" HighMass="%f" SamplesPerAMU="%f" Units="%[^"]" PiraniPressure="%f" TotalPressure="%f" FilamentStatus="%f" Sample="%f">','headerLines',1,'expChars','e');
fseek(fileid, 0, 'bof');
data=textscan(fileid,'<Sample Value="%f"/>','expChars','e','headerLines',2);
if isnumeric(data{1,1})==1
    dataout(:,1)=(0:(length(data{1,1})-1)).'/(length(data{1,1}))*(header{2}-header{1})+header{1};
    dataout(:,2)=data{1,1};
end
varargout{1}=header;
fclose(fileid);
end
function varargout = tracehit(varargin)

if nargin==0
    [FileName,PathName] = uigetfile({'*.*','All Files (*.*)'},'Choose a trace.');
    trace=load(strcat(PathName,FileName));
    if exist('TOF_bins.txt')==2
        bin_file='TOF_bins.txt';
    else [FileName1,PathName1] = uigetfile({'Text Files (*.txt)'},'Choose a bin file.');
        bin_file=[PathName1 FileName1];
    end
    answer = inputdlg('Threshhold?');
    thrs=str2num(answer{1});
elseif nargin==1
    if ischar(varargin{1})
        trace=load(varargin{1});
    elseif isnumeric(varargin{1})
        trace=varargin{1};    
    else 'Wrong type of argument.'
        return;
    end
    if exist('TOF_bins.txt')==2
        bin_file='TOF_bins.txt';
    else [FileName1,PathName1] = uigetfile({'Text Files (*.txt)'},'Choose a bin file.');
        bin_file=[PathName1 FileName1];
    end
    thrs=150;
elseif nargin==2
    if ischar(varargin{1})
        trace=load(varargin{1});
    elseif isnumeric(varargin{1})
        trace=varargin{1};    
    else 'Wrong type of argument.'
        return;
    end
    if exist('TOF_bins.txt')==2
        bin_file='TOF_bins.txt';
    else [FileName1,PathName1] = uigetfile({'Text Files (*.txt)'},'Choose a bin file.');
        bin_file=[PathName1 FileName1];
    end
    thrs=varargin{2};
elseif nargin==3
    if ischar(varargin{1})
        trace=load(varargin{1});
    elseif isnumeric(varargin{1})
        trace=varargin{1};    
    else 'Wrong type of argument.'
        return;
    end
    bin_file=varargin{3};
    thrs=varargin{2};
end

% tic;
fid=fopen(bin_file);
temp = textscan(fid, '%s %d %d %d %d');
fclose(fid);
% t=toc
for ind1=1:length(temp{1})
    peaknamematrix{ind1}=temp{1}{ind1};
    S(ind1).name=temp{1}{ind1};
    peakmatrix(ind1,:)=[temp{2}(ind1) temp{3}(ind1)];
    backgroundmatrix(ind1,:)=[temp{4}(ind1) temp{5}(ind1)];
end
N=length(trace);
% cut=[2:150 465:590 642:868 1847:2113 (ceil(N/2)+1):N];
a2=trace-mean(trace(end-100:end));
% aft=fft(a2);
% aft2=aft;
% aft2(cut)=0;
% aft2=2*aft2;
% aift=abs(ifft(aft2));
if max(a2(peakmatrix(1,1):peakmatrix(1,2)))>=thrs % do we have an electron hit on the trace?
    counts(1)=1;
    for ind2=2:size(peakmatrix,1) % if yes, check the ion bins too
        if max(a2(peakmatrix(ind2,1):peakmatrix(ind2,2)))>=thrs; % just check if the signal goes above threshold
            counts(ind2)=1;
        else
            counts(ind2)=0;
        end
    end
    detected=unique((1:length(peaknamematrix)).*counts);
    detected(detected==0)=[];
    varargout{1}=detected;
else
    varargout{1}=[];
end
end
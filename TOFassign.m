function TOFassign(varargin)

if nargin==1
%     [PathName, name, ext, versn] = fileparts(varargin{1});
%     PathName=strcat(PathName,'\');
%     FileName=strcat(name,ext);
    if ischar(varargin{1})==1
        if isempty(strfind(varargin{1},'.bin'))==0
            S=loadbinarytrace(varargin{1});
        else S=importdata(varargin{1});
        end
    elseif isfloat(varargin{1})==1
        S=varargin{1};
    end
elseif nargin==0
    [FileName,PathName] = uigetfile({'*.dat*','Trace file';'*.txt','Text File';'*.bin','Binary File';'*.*','All Files (*.*)'},'Pick a trace file');
    if isempty(strfind(FileName,'.bin'))==0
        S=loadbinarytrace(strcat(PathName,FileName));
    else S=importdata(strcat(PathName,FileName));
    end
end
% S=importdata(strcat(PathName,FileName));
% S=uiimport('-file',{'*.scope','OceanOptics spectrum file';'*.txt','text data';'*.*','All Files (*.*)'},'Pick a spectrum file');
if isstruct(S)==1
    data=map2colvec(S.data);
elseif isnumeric(S)==1
    data=map2colvec(S);
end
x0=map2colvec(1:length(data));
set_mass_peak_binsP;
h1=figure;h2=plot(data,'b.-');
% h1=figure;h2=plot(TOFcalib(x0,A,t0),data,'b.-');
for ind1=1:length(peaknamematrix)
    x=round(sum(x0(peakmatrix(ind1,1):peakmatrix(ind1,2)).*data(peakmatrix(ind1,1):peakmatrix(ind1,2)))/sum(data(peakmatrix(ind1,1):peakmatrix(ind1,2))));
    y=data(x);
    text(x,y,peaknamematrix{ind1},'fontsize',14)
    xlabel('TOF delay [arb. units]')
    ylabel('Fragment yields [arb. units]')
    setfigP
%     ax2=axes('position',get(gca,'position'),'Xaxislocation','top','Yaxislocation','right','Color','none','Ycolor',[.8,.8,.8],'Ytick',[])
%     set(get(ax2,'children'),'xdata',TOFcalib(x0,A,t0))
%     xlabel(ax2,'Mass [atomic units]','fontsize',14)    
end
saveas(h2,'peak_assignments','fig');
saveas(h2,'peak_assignments','png');
close(h1)
end
function textout = cellcat(varargin)

if nargin==0
    [FileName,PathName] = uigetfile({'*.txt','Text File';'*.*','All Files (*.*)'},'Pick a textfile');
    fileID=fopen(strcat(PathName,FileName));
    textdata=textscan(fileID,'%s');
elseif nargin>0
    textdata=varargin{1};
end
textout=[];
if iscell(textdata)
    for ind1=1:length(textdata{1})
        textout=[textout textdata{1}{ind1}];
        textout=[textout ' '];
    end
elseif ischar(textdata)
    textout=textdata;
else 'No valid textdata found.'
end
end
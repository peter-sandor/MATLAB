function out=getdispname(varargin)
if nargin==0
    out=str2num(get(gco,'DisplayName'));
else
    out=str2num(cell2mat(get(varargin{1},'DisplayName')));
end
end
function varargout = dfigure(varargin)
if nargin==0
    hfig = figure('windowstyle','docked');
else hfig = figure('windowstyle','docked',varargin{1},varargin{2});
end
varargout{1}=hfig;
end

function varargout = getdata(varargin)

if nargin==0
    hndl_ax=gca;
elseif nargin==1;
    hndl_ax=varargin{1};
end
x0=get(get(hndl_ax,'Children'),'xdata');
y0=get(get(hndl_ax,'Children'),'ydata');
Nx=length(x0);
if iscell(x0)
% The order in which the curves are plotted is reversed.
% Fix it below.
    for ind1=1:Nx
        x1{ind1}=x0{Nx-ind1+1};
        y1{ind1}=y0{Nx-ind1+1};
    end
    varargout{1}=x1;
    varargout{2}=y1;
else
    varargout{1}=x0;
    varargout{2}=y0;
end
if isprop(get(hndl_ax,'Children'),'cdata')
    varargout{3}=get(get(hndl_ax,'Children'),'cdata');
elseif isprop(get(hndl_ax,'Children'),'zdata')
    varargout{3}=get(get(hndl_ax,'Children'),'zdata');
else
    varargout{3}=[];
end
if isempty(ishandle(legend))==0
    strct_legend=legend;
    cell_legend0=strct_legend.String;
%     for ind1=1:Nx
%         cell_legend1{ind1}=cell_legend0{Nx-ind1+1};
%     end
    if isempty(varargout{3})==0
        varargout{3}=cell_legend0;
    else varargout{4}=cell_legend0;
    end
end
end
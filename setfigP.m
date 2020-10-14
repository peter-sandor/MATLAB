function setfigP(varargin)
% This script sets all property values of the selected axes to 
% optimal values for papers and (small) prints.

% written by Florian 05/08/04
% edited by Peter 2012 spring
% optional input: axes handle
if nargin==0
    h=gca;
    fntsz=16;
    lnthck=2;
elseif nargin==1
    if ishandle(varargin{1})
        h=varargin{1};
        fntsz=16;
        lnthck=2;
    else
        h=gca;
        fntsz=varargin{1}(1);
        lnthck=varargin{1}(2);
    end
elseif nargin==2
    h=varargin{1};
    fntsz=varargin{2}(1);
    lnthck=varargin{2}(2);
end
set(h,'fontname','arial');
set(h,'fontsize',fntsz);
% set(h,'fontweight','bold');
set(h,'gridlinestyle','--');
set(h,'linewidth',lnthck);
set(get(h,'xlabel'),'fontsize',fntsz)
set(get(h,'ylabel'),'fontsize',fntsz)
set(get(h,'title'),'fontsize',fntsz)
if isempty(get(h,'children'))==0
    if iscell(get(h,'children'))==1
        set(cell2mat(get(h,'children')),'linewidth',lnthck);
    else
        if ~isprop(get(h,'children'),'CDataMapping')
            set(get(h,'children'),'linewidth',lnthck);
        end
    end
end
if isempty(findobj(gcf,'Type','axes','Tag','legend'))==0
    set(findobj(gcf,'Type','axes','Tag','legend'),'fontsize',fntsz)
end
end


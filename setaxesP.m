function setaxesP(varargin)
% This script sets all property values of the selected axes to 
% optimal values for papers and (small) prints.

% written by Florian 05/08/04
% edited by Peter 2012 spring
if nargin==0
    h=gca;
elseif nargin>0
    h=varargin{1};
end
set(h,'fontname','Times');
set(h,'fontsize',14);
% set(h,'fontweight','bold');
set(h,'gridlinestyle','--');
set(h,'linewidth',3);
if isempty(get(h,'children'))==0
    if iscell(get(h,'children'))==1
        hndl_lines=get(h,'children');
        for ind1=1:length(hndl_lines)
            set(hndl_lines{ind1},'linewidth',2);
        end
    else set(get(h,'children'),'linewidth',2);
    end
end
if isempty(findobj(gcf,'Type','axes','Tag','legend'))==0
    set(findobj(gcf,'Type','axes','Tag','legend'),'fontsize',12)
end
end


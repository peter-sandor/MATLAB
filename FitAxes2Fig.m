function FitAxes2Fig(varargin)
if nargin==0
    hfig=gcf;
    hax=findobj(hfig,'type','axes');
else
    hax=varargin{1};
end

% Default settings for no xlabel, no ylabel, to title,
% xaxislocation=bottom, yaxislocation=left.
% xstart=0.07;
% xend=0.98;
% ystart=0.07;
% yend=0.98;

struct_xlabel=get(hax,'xlabel');
struct_ylabel=get(hax,'ylabel');
struct_title=get(hax(1),'title');

if length(hax)==1
    if strcmp(get(hax,'xaxislocation'),'top')
        ystart=0.02;
        if isempty(struct_xlabel.String)
            yend=0.93;
        else 
            yend=0.86;
        end
    else
        yend=0.98;
        if isempty(struct_xlabel.String)
            ystart=0.07;
        else 
            ystart=0.14;
        end
    end
    
    if strcmp(get(hax(1),'yaxislocation'),'right')
    xstart=0.02;
        if isempty(struct_ylabel.String)
            xend=0.93;
        else 
            xend=0.88;
        end
    else
    xend=0.98;
        if isempty(struct_ylabel.String)
            xstart=0.07;
        else 
            xstart=0.12;
        end
    end

else
    if strcmp(get(hax(1),'yaxislocation'),'right')
    xstart=0.02;
        if isempty(struct_ylabel{1}.String)
            xend=0.93;
        else 
            xend=0.88;
        end
    else
    xend=0.98;
        if isempty(struct_ylabel{1}.String)
            xstart=0.07;
        else 
            xstart=0.12;
        end
    end
end

if ~isempty(struct_title.String)
    yend=yend-0.04;
end

axpos=get(hax,'position');
if length(hax)==1
    axpos2=[xstart ystart xend-xstart yend-ystart];
    set(hax,'position',axpos2);
else
    for ind1=1:length(hax)
        axpos2{ind1}=[xstart axpos{ind1}(2) xend-xstart axpos{ind1}(4)];
        set(hax(ind1),'position',axpos2{ind1});
    end
end
end
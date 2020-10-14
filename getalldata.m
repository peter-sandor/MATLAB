function data = getalldata(varargin)

if nargin==0
    hfig=gcf;
elseif nargin==1;
    hfig=varargin{1};
end

hAllAxes = findobj(hfig,'type','axes');
hLeg = findobj(hfig,'type','legend');
% hAxes = setdiff(hAllAxes,hLeg); % All axes which are not legends

for ind1=1:length(hAllAxes)
    hndl_child=get(hAllAxes(ind1),'Children');
    hndl_line=findobj(hndl_child,'type','line');
    data(ind1).xdata=get(hndl_line,'xdata');
    data(ind1).ydata=get(hndl_line,'ydata');
    if isprop(hndl_child,'cdata')
        data(ind1).zdata=get(hndl_child,'cdata');
        data(ind1).xdata=get(hndl_child,'xdata');
        data(ind1).ydata=get(hndl_child,'ydata');
    elseif isprop(hndl_child,'zdata')
        data(ind1).zdata=get(hndl_child,'zdata');
        data(ind1).xdata=get(hndl_child,'xdata');
        data(ind1).ydata=get(hndl_child,'ydata');
    end
    if isempty(ishandle(hLeg))==0
        str_legend=legend;
        data(ind1).legend=map2colvec(flipdim(str_legend.String,2));
    end
end
end
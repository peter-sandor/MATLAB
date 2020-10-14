function varargout = getplot(varargin)

if nargin==0
    hndl1=gca;
elseif nargin==1
    hndl1=varargin{1};
end

[struct_plot.x,struct_plot.y,struct_plot.z]=getdata(hndl1);

if isempty(ishandle(get(hndl1,'Title')))==0
    struct_plot.title=get(get(hndl1,'Title'),'String');
end
if isempty(ishandle(get(hndl1,'xlabel')))==0
    struct_plot.xlabel=get(get(hndl1,'xlabel'),'String');
end
if isempty(ishandle(get(hndl1,'ylabel')))==0
    struct_plot.ylabel=get(get(hndl1,'ylabel'),'String');
end
    
% hndl4=findall(hndl1);
% for ind1=2:length(hndl4)
%     set(hndl4(ind1),'parent',hax1)
% end

struct_plot.xlim=get(hndl1,'xlim');
struct_plot.ylim=get(hndl1,'ylim');
struct_plot.clim=get(hndl1,'clim');
struct_plot.xscale=get(hndl1,'xscale');
struct_plot.yscale=get(hndl1,'yscale');
varargout{1}=struct_plot;
end
function varargout = getdata

varargout{1}=get(get(gca,'Children'),'xdata');
varargout{2}=get(get(gca,'Children'),'ydata');
if isprop(get(gca,'Children'),'cdata')
    varargout{3}=get(get(gca,'Children'),'cdata');
elseif isprop(get(gca,'Children'),'zdata')
    varargout{3}=get(get(gca,'Children'),'zdata');
end
end
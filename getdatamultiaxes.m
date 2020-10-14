function varargout = getdatamultiaxes

h=get(gcf,'children'); % get handles for all axes belonging to current figure
h2=[];
for ind1=1:length(h)
    h2 = [h2; findobj(h(ind1),'Type','line')];
end
varargout{1}=get(h2,'xdata');
varargout{2}=get(h2,'ydata');
% for ind2=1:length(h2)
%     if isprop(h2(ind2),'cdata')
%         varargout{3}=get(h2(ind2),'cdata');
%     elseif isprop(h2(ind2),'zdata')
%         varargout{3}=get(h2(ind2),'zdata');
%     end
% end
if isempty(ishandle(findobj(gcf,'Type','axes','Tag','legend')))==0
    if isempty(varargout{3})==0
        varargout{3}=map2colvec(flipdim(get(findobj(gcf,'Type','axes','Tag','legend'),'String'),2));
    else varargout{4}=map2colvec(flipdim(get(findobj(gcf,'Type','axes','Tag','legend'),'String'),2));
end
end
function replace_axes(x,y)
% This code replaces the x-axis data on the current figure ('origx') with the one supplied in the argument ('x')
% There is a built-in check to see if the vector lengths match, otherwise the figure would get screwed up.

if ~isempty(x)
    origx=get(get(gca,'Children'),'xdata');
    if iscell(origx)==1
        origx=origx{1};
    end
    if length(origx)==length(x)
        set(get(gca,'Children'),'xdata',x);
    else disp('Make sure the x-axis has the proper number of points!');
    end
end

if ~isempty(y)
    origy=get(get(gca,'Children'),'ydata');
    if iscell(origy)==1
        origy=origy{1};
    end
    if length(origy)==length(y)
        set(get(gca,'Children'),'ydata',y);
    else disp('Make sure the y-axis has the proper number of points!');
    end
end
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
end
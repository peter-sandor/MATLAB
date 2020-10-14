function [x,y,z] = getdata3d

x=get(get(gca,'Children'),'xdata');
y=get(get(gca,'Children'),'ydata');
z=get(get(gca,'Children'),'cdata');

end
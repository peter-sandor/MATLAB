function x0 = ip_point(x,y,y0)

% finds the abscissa point x0 for the value of y0, obtained using linear
% interpolation between x1 --> y1 and x2 --> y2
% latter ones are stored in vectors x and y
% it seems data type 'double' stores 15 digits
if (y0<=max(y))&&(y0>=min(y))
    y=y-y0;
    x0 = x(1) + abs(y(1))/(abs(y(1))+abs(y(2)))*(x(2)-x(1));
else 'y0 is not not within limits'
end
end
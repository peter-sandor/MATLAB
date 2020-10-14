function dydx = num_diff(x,y)

N = length(x);
dydx0 = diff(y)./diff(x);
for ind1 = 1:N-1
    x0(ind1) = (x(ind1+1) + x(ind1))/2;
end
dydx = interp1(x0,dydx0,x(1:end),'linear','extrap');

% dydx = interp1(x0,dydx0,x(2:end-1),'linear');
% dydx = map2colvec([(dydx(2)-dydx(1))/(x0(2)-x0(1))*(x(1)-(x0(2)+x0(1))/2) dydx dydx(end)]);
end
function pos_out = genpos_for_knifeedge(pos_center,pos_span,pos_sigma)
% input arguments in [um]
N = 100;
min_stepsize = 1;
max_stepsize = 10;
vec_base = pos_center-pos_span/2:pos_span/(N-1):pos_center+pos_span/2;
pos_density = round(max_stepsize - (max_stepsize-min_stepsize)*exp(-(vec_base-pos_center).^2/pos_sigma^2));
pos_out(1) = round(pos_center-pos_span/2);
stepsize(1) = 10;
ind1 = 2;
while pos_out(ind1-1)<max(vec_base)
    stepsize(ind1) = round(interp1(vec_base,pos_density,pos_out(ind1-1)));
    pos_out(ind1) = pos_out(ind1-1) + stepsize(ind1);
    ind1 = ind1+1;
end
pos_out = map2colvec(pos_out*1000);
end
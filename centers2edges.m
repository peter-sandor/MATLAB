function x_edge = centers2edges(x_center)
% This function generates histogram bin edges from bin centers.
% For use with 'histcounts' and 'bar' functions.

x_center = map2colvec(x_center);
N = length(x_center);
x_diff = diff(x_center);
x_edge = x_center(1:N-1) + x_diff/2;
x_edge = [x_center(1)-x_diff(1)/2; x_edge; x_center(N)+x_diff(N-1)/2];
end
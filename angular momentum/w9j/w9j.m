%-------------------------------------------------------------------------%
%Filename:  w9j.m
%Author:    Oliver Johnson
%Date:      6/7/2011
%
% W9J computes the Wigner 9-j coefficients.
%
% Inputs:
%   a,b,c,d,e,f,g,h,j - Scalar arguments to the 9-j symbol from the form:
%
%                                   {a b c}
%                                   {d e f}
%                                   {g h j}
%
% Outputs:
%   W - A scalar giving the required Wigner 9-j coefficient.
%
% Note: Tested on first 87 cases from Table 10.13 of Ref. [1].
%
% [1] Varshoalovich; Quantum Theory of Angular Momentum. (1988). p. 340.
% [2] Format of xlims taken from:
%     http://www.mathworks.com/matlabcentral/fileexchange/20619
% [3] Weisstein, Eric W. "Triangular Inequalities."
%     From MathWorld--A Wolfram Web Resource.
%     http://mathworld.wolfram.com/TriangularInequalities.html
%-------------------------------------------------------------------------%

function W = w9j(a,b,c,d,e,f,g,h,j)

assert(isscalar(a) && isscalar(b) && isscalar(c) && ...
    isscalar(d) && isscalar(e) && isscalar(f) && ...
    isscalar(g) && isscalar(h) && isscalar(j),'All inputs must be scalars.')

%---check conditions---%
args = [a b c d e f g h j];
if any(args < 0) || any(mod(2*args,1) ~= 0)
    
    error('All arguments to Wigner 9-j symbol must be integer or half-integer non-negative numbers.')
    
elseif ~triangular_cond(a,b,c) || ~triangular_cond(d,e,f) || ...
        ~triangular_cond(g,h,j) || ~triangular_cond(a,d,g) || ...
        ~triangular_cond(b,e,h) || ~triangular_cond(c,f,j)
    
    W = 0;
    return
    
end

%---compute valid x values for summation---%
xlims(1) = a+j;
xlims(2) = b+f;
xlims(3) = d+h;
xlims(4) = a-j;
xlims(5) = b-f;
xlims(6) = d-h;

%---compute W---%
W = 0;
for x = max(xlims(4:6)):min(xlims(1:3))
    
    W = W + (-1)^(2*x)*(2*x+1)*...
        w6j(a,b,c,f,j,x)*w6j(d,e,f,b,x,h)*w6j(g,h,j,x,a,d);
    
end

end

function tf = triangular_cond(a,b,c)

tf = (c >= abs(a-b)) & (c <= a+b);

end
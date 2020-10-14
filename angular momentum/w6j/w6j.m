%-------------------------------------------------------------------------%
%Filename:  w6j.m
%Author:    Oliver Johnson
%Date:      6/7/2011
%
% W6J computes the Wigner 6-j coefficients.
%
% Inputs:
%   a,b,c,d,e,f - Scalar arguments to the 6-j symbol from the form:
%
%                                   {a b c}
%                                   {d e f}
%
% Outputs:
%   W - A scalar giving the required Wigner 6-j coefficient.
%
% Note: Tested on all 2264 cases given at [3], and results were within eps.
%
% [1] Varshoalovich; Quantum Theory of Angular Momentum. (1988). p. 293.
% [2] Messiah; Quantum Mechanics, Vol. 2. (1962). p.1063.
% [3] http://www.strw.leidenuniv.nl/~mathar/progs/6jSymb
% [4] Format of nlims taken from:
%     http://www.mathworks.com/matlabcentral/fileexchange/20619
%-------------------------------------------------------------------------%

function W = w6j(a,b,c,d,e,f)

assert(isscalar(a) && isscalar(b) && isscalar(c) && isscalar(d) && isscalar(e) && isscalar(f),'All inputs must be scalars.')

%---check conditions---%
if ~triangular_cond(a,b,c) || ~triangular_cond(c,d,e) || ...
        ~triangular_cond(a,e,f) || ~triangular_cond(b,d,f) || ...
        mod(a+b+c,1) ~= 0 || mod(c+d+e,1) ~= 0 || ...
        mod(a+e+f,1) ~= 0 || mod(b+d+f,1) ~= 0
    
    W = 0;
    return
    
end

%---compute valid n values for summation---%
nlims(1) = 0;
nlims(2) = a+b+c;
nlims(3) = c+d+e;
nlims(4) = a+e+f;
nlims(5) = b+d+f;
nlims(6) = a+b+d+e;
nlims(7) = a+c+d+f;
nlims(8) = b+c+e+f;
n = max(nlims(1:5)):min(nlims(6:8));

%---check for stability---%
if a+b-c > 21 || a-b+c > 21 || -a+b+c > 21 || a+b+c+1 > 21 || ...
        c+d-e > 21 || c-d+e > 21 || -c+d+e > 21 || c+d+e+1 > 21 || ...
        a+e-f > 21 || a-e+f > 21 || -a+e+f > 21 || a+e+f+1 > 21 || ...
        b+d-f > 21 || b-d+f > 21 || -b+d+f > 21 || b+d+f+1 > 21 || ...
        any(n+1 > 21) || ...
        any(n-nlims(2) > 21) || any(n-nlims(3) > 21) || ...
        any(n-nlims(4) > 21) || any(n-nlims(5) > 21) || ...
        any(nlims(6)-n > 21) || any(nlims(7)-n > 21) || ...
        any(nlims(8)-n > 21)
    
    warning('The argument to one or more of the factorials used in the computation of the requested Clebsch-Gordan coefficient is greater than 21, this can result in inaccuracies (see Matlab documentation for the FACTORIAL function).') %#ok<WNTAG>
    
end

%---compute coefficient---%
W = del(a,b,c)*del(c,d,e)*del(a,e,f)*del(b,d,f)*...
    sum((((-1).^n).*factorial(n+1))./...
    (factorial(n-nlims(2)).*factorial(n-nlims(3)).*...
    factorial(n-nlims(4)).*factorial(n-nlims(5)).*...
    factorial(nlims(6)-n).*factorial(nlims(7)-n).*...
    factorial(nlims(8)-n)));

end

function tri = del(a,b,c)

tri = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1)); %Eq. 8.2(1)

end

function tf = triangular_cond(a,b,c)

tf = (c >= abs(a-b)) & (c <= a+b);

end
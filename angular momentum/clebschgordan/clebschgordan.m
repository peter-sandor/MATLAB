%-------------------------------------------------------------------------%
%Filename:  clebschgordan.m
%Author:    Oliver Johnson
%Date:      6/7/2011
%
% CLEBSCHGORDAN computes Clebsch-Gordan coefficients.
%
% Inputs:
%   j1 - A scalar giving the first total angular momentum.
%   m1 - A scalar giving the projection of the first total angular
%        momentum.
%   j2 - A scalar giving the second total angular momentum.
%   m2 - A scalar giving the projection of the second total angular
%        momentum.
%   j  - A scalar giving the coupled total angular momentum.
%   m  - A scalar giving the projection of the coupled total angular
%        momentum.
%
% Outputs:
%   C - A scalar giving the required Clebsch-Gordan coefficient.
%
% Note: Tested on all 7392 cases given at [3], and results were within eps.
%
% [1] G. Racah; Phys. Rev., 1942, vol. 62, pp. 438-462.
% [2] Varshoalovich; Quantum Theory of Angular Momentum. (1988). Chapter 8.
% [3] http://www.strw.leidenuniv.nl/~mathar/progs/CGord
%-------------------------------------------------------------------------%

function C = clebschgordan(j1,m1,j2,m2,j,m)

assert(isscalar(j1) && isscalar(m1) && isscalar(j2) && isscalar(m2) && isscalar(j) && isscalar(m),'All inputs must be scalars.')

%---check conditions---%
if j1 < 0 || j2 < 0 || j < 0 || ...
        mod(2*j1,1) ~= 0 || mod(2*j2,1) ~= 0 || mod(2*j,1) ~= 0 || ...
        mod(2*m1,1) ~= 0 || mod(2*m2,1) ~= 0 || mod(2*m,1) ~= 0 || ...
        abs(m1) > j1 || abs(m2) > j2 || abs(m) > j || ...
        j1+m1 < 0 || j2+m2 < 0 || j+m < 0 || j1+j2+j < 0 ||...
        mod(j1+m1,1) ~= 0 || mod(j2+m2,1) ~= 0 || mod(j+m,1) ~= 0 || ...
        mod(j1+j2+j,1) ~= 0
    
    error(sprintf('Clebsch-Gordan coefficient only defined if: \n 1. j1, j2, j are integer or half-integer non-negative numbers. \n 2. m1, m2, m are integer or half-integer numbers. \n 3. abs(m1)<=j1, abs(m2)<=j2, abs(m)<=j \n 4. j1+m1, j2+m2, j+m, j1+j2+j are integer non-negative numbers.')) %#ok<SPERR>
   
elseif m1+m2-m ~= 0 || j < abs(j1-j2) || j > j1+j2
    
    C = 0;
    return
    
end

%---compute valid k values for summation---%
k = max([0,j2-j-m1,j1-j+m2]):min([j1+j2-j,j1-m1,j2+m2]);

%---check for stability---%
if j+j1-j2 > 21 || j+j2-j1 > 21 || j1+j2-j > 21 || j1+j2+j+1 > 21 || ...
        j+m > 21 || j-m > 21 || j1+m1 > 21 || j1-m1 > 21 || ...
        j2+m2 > 21 || j2-m2 > 21 || any(k > 21) || ...
        any(j1+j2-j-k > 21) || any(j1-m1-k > 21) || ...
        any(j2+m2-k > 21) || any(j-j2+m1+k > 21) || any(j-j1-m2+k > 21)
    
%     warning('The argument to one or more of the factorials used in the computation of the requested Clebsch-Gordan coefficient is greater than 21, this can result in inaccuracies (see Matlab documentation for the FACTORIAL function).') %#ok<WNTAG>
    
end

%---compute coefficient---%
C = sqrt((2*j+1)*factorial(j+j1-j2)*factorial(j+j2-j1)*factorial(j1+j2-j)/factorial(j1+j2+j+1))*...
    sqrt(factorial(j+m)*factorial(j-m)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2))*...
    sum(((-1).^k)./(factorial(k).*factorial(j1+j2-j-k).*factorial(j1-m1-k).*factorial(j2+m2-k).*factorial(j-j2+m1+k).*factorial(j-j1-m2+k)));
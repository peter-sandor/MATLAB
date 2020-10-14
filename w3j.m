%-------------------------------------------------------------------------%
%Filename:  w3j.m
%Author:    Oliver Johnson
%Date:      6/7/2011
%
% W3J computes the Wigner 3-j coefficients.
%
% Inputs:
%   j1 - A scalar giving the first total angular momentum.
%   m1 - A scalar giving the projection of the first total angular
%        momentum.
%   j2 - A scalar giving the second total angular momentum.
%   m2 - A scalar giving the projection of the second total angular
%        momentum.
%   j3 - A scalar giving the third total angular momentum.
%   m3 - A scalar giving the projection of the third total angular
%        momentum.
%
% Outputs:
%   W - A scalar giving the required Wigner 3-j coefficient.
%
% [1] http://en.wikipedia.org/wiki/3-jm_symbol.
%-------------------------------------------------------------------------%

function W = w3j(j1,j2,j3,m1,m2,m3)

W = (-1)^(j1-j2-m3)*clebschgordan(j1,m1,j2,m2,j3,-m3)/sqrt(2*j3+1);
function [coeffs,LegPoly] = expand_Legendre(x,f,N)
% This function expands f(x), on range -1<=x<=1 on the basis of the first N
% Legendre-polynomials, and returns the expansion coefficients and the polynomials.
x=map2colvec(x);
f=map2colvec(f);
dx=[diff(x); (max(x)-min(x))/length(x)];
for ind1=1:N
    temp = legendre(ind1-1,x).';
    LegPoly(:,ind1)=temp(:,1);
    coeffs(ind1)=(2*(ind1-1)+1)/2*sum(LegPoly(:,ind1).*dx.*f,1);
end
coeffs=map2colvec(coeffs);
end
function out = wigner_d(theta,JMK)
%Calculates Wigner_D functions, the eigenfunctions for rotational states
%of symmetric top molecules (only the theta part of the eigenfunction is
%computed, the complex exponentials which comprise the phi and chi
%contributions are omitted)
%
%Enter quantum numbers for eigenstate of interest
J=JMK(1);
K=JMK(3);
M=JMK(2);
%Enter number of theta points (between 0 and pi)
% numpoints=1000;
% thetastep=pi/(numpoints+1);
%compute d
% theta=(0:thetastep:pi);
c=cos(theta/2);
s=sin(theta/2);
out=0*s;
prefactor=sqrt(J+0.5)*((-1)^(J-M))*sqrt(factorial(J+M)*factorial(J-M)*factorial(J+K)*factorial(J-K));
check1=[(J-M),(J-K)];
maxN=min(check1);
check2=[(K+M), 0];
minN=-min(check2);
for N=minN:maxN
   out=out + ((-1)^N)*(c.^(2*N+M+K)).*(s.^(2*J-2*N-M-K))/(factorial(K+M+N)*factorial(J-M-N)*factorial(J-K-N)*factorial(N));
end
out=out*prefactor;
end


function func_d = wigner_d_P(theta,maxJ)

Ntheta=length(theta);
func_d = zeros([Ntheta maxJ+1 2*maxJ+1 2*maxJ+1]);
Jvals=0:maxJ;
cos_theta_d2=map2colvec(cos(theta/2));
sin_theta_d2=map2colvec(sin(theta/2));
for indJ=1:maxJ+1
    J=Jvals(indJ);
    Mvals=-J:J;
    offset=maxJ+1-indJ;
    for indm1=1:(2*J+1)
        M=Mvals(indm1);
        for indm2=1:(2*J+1)
            K=Mvals(indm2);
%             prefactor=sqrt(factorial(J+M)*factorial(J-M)/factorial(J+K)/factorial(J-K));
            prefactor=(-1)^(J-M)*sqrt((J+1/2)*factorial(J+M)*factorial(J-M)/factorial(J+K)/factorial(J-K));
            max_s=min([J-K J-M]);
            min_s=max([-K-M 0]);
%             min_s=-min([K+M 0]);
            for s=min_s:max_s
                func_d(:,indJ,indm1+offset,indm2+offset) = func_d(:,indJ,indm1+offset,indm2+offset) + prefactor*(-1)^s*factorial(J-K)*factorial(J+K)/factorial(J-K-s)/factorial(J-M-s)/factorial(s)/factorial(M+K+s)*cos_theta_d2.^(2*s+K+M).*sin_theta_d2.^(2*J-2*s-M-K);
            end
        end
    end
end
end
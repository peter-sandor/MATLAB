Ndelay=100;
maxJ=50;
Ntheta=200;
theta=0:pi/(Ntheta-1):pi;
legendre_poly=LegendreP(cos(theta),maxJ);
Ylm_reduced=zeros(size(legendre_poly));
ang_distr=zeros([Ntheta Ndelay]);
for indM=1:54%2*maxJ+1;
    if indM==53
        1;
    end
    for indJ=1:maxJ+1
            Jval=indJ-1;
            Mval=indM-(maxJ+1);
            if Jval>=abs(Mval)
                factorM=1/prod([Jval-Mval+1:1:Jval+Mval]);
                Ylm_reduced(:,indJ,indM)=(-1)^Mval*sqrt((2*Jval+1)/4/pi*factorM)*legendre_poly(:,indJ,indM);
            end
    end
end
%%
% load simresult;
phi_step=2*pi/99;
phi=0:2*pi/(99):2*pi;
theta_ext=permute(extend(simresult.theta,length(simresult.delay)),[1 2]);
indJ=30;
indM=52;
Ylm=sph_harmonics(simresult.theta,phi,indJ-1,indM-51);
sum(sum(Ylm.*cos(theta_ext).^2.*sin(theta_ext).*conj(Ylm)))*thstep*2*pi/(99)

sum(conj(simresult.Ylm_red(:,indJ,indM)).*simresult.Ylm_red(:,indJ,indM).*cos(map2colvec(simresult.theta)).^2.*sin(map2colvec(simresult.theta)))*thstep*2*pi
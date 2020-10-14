function [Rbins polarhist] = polar_binner(theta,rho,Nrad)

% rho: Ntheta x N array
Ntheta=length(theta);
Rmax=1.05*max(max(rho));
Ncurve=size(rho,2);
samplingfactor=Nrad/Rmax;
Rstep=Rmax/(Nrad-1);
Rbins=map2colvec(0:Rstep:Rmax);
polarhist=zeros([Ntheta Nrad]);

for ind1=1:Ncurve
    for ind2=1:Ntheta
        rhoind=sum(double(rho(ind2,ind1)>Rbins));
        if rhoind==0
            rhoind=1;
        end
        polarhist(ind2,rhoind)=polarhist(ind2,rhoind)+1;
    end
end
Rbins=Rbins+Rstep/2;
end
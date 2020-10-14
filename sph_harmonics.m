function out = sph_harmonics(theta,phi,l,m)

if abs(m)>l || mod(l,1)~=0 || mod(m,1)~=0
    out=[];
    disp('l and m should be integers and |m|<=l');
else
    Ntheta=length(theta);
    Nphi=length(phi);
    for ind1=1:Ntheta
        temp(:,:,ind1)=legendre(l,cos(theta(ind1)));
        if m<0
            lgndr(ind1)=(-1)^abs(m)*factorial(l-abs(m))/factorial(l+abs(m))*map2colvec(squeeze(temp(abs(m)+1,:,ind1)));
        else
            lgndr(ind1)=map2colvec(squeeze(temp(m+1,:,ind1)));
        end
    end
    out=(-1)^m*sqrt((2*l+1)/4/pi*factorial(l-m)/factorial(l+m))*map2colvec(lgndr)*exp(1i*m*map2rowvec(phi));
end

end
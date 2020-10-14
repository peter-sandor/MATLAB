function [nrhits,coords] = getshots(imgin,thrs)

N=size(imgin);
% N=flipdim(N,2); % convert to [x y] = [column row] format
% chop up the image to subimages no smaller than 60x60 but no larger than
% 100x100.
Nfac1=factor(N(1));
Nfac2=factor(N(2));
cpr1=cumprod(Nfac1);
cpr2=cumprod(Nfac2);
M=[min(cpr1(cpr1>=60)) min(cpr2(cpr2>=60))];
nrchunk=N./M;
nrhits=0;
coords=[];
for ind1=0:nrchunk(1)-1
    for ind2=0:nrchunk(2)-1
        [a,b]=centroid(imgin(ind1*M(1)+1:(ind1+1)*M(1),ind2*M(2)+1:(ind2+1)*M(2)));
        b=b+extend([ind2*M(2)+1:(ind2+1)*M(2) ind1*M(1)+1:(ind1+1)*M(1)],size(b,1));
        nrhits=nrhits+a;
        coords=cat(1,coords,b);
    end
end
end
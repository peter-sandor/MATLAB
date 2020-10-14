function out = meanw(A,w,dim)
% This code averages the contents array 'A' along dimension 'dim' using
% weights 'w'. E.g. if 'A' is an 3x4x5 array, and 'w' is a 1x4 array then
% the call 'meanw(A,w,2)' will result in the 3x5 array 'out'.

Np=size(A);
if Np(dim)==length(w);
    cpr=cumprod(Np([1:dim-1 dim+1:length(Np)]));
    cpr2=cumprod(Np);
    out=zeros([cpr(end) 1]);
    wsum=sum(w);
    for ind1=1:cpr2(end)
        indw=ind2subP(Np,ind1);
        out(sub2indP(Np([1:dim-1 dim+1:end]),indw([1:dim-1 dim+1:end])))=out(sub2indP(Np([1:dim-1 dim+1:end]),indw([1:dim-1 dim+1:end])))+A(ind1)*w(indw(dim))/wsum;
    end
    out=reshape(out,Np([1:dim-1 dim+1:length(Np)]));
else 'Subscripted dimensions dont match.'
end
end
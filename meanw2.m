function out = meanw2(A,w,dim)
% This code averages the contents array 'A' along dimension 'dim' using
% weights 'w'. E.g. if 'A' is an 3x4x5 array, and 'w' is a 1x4 array then
% the call 'meanw(A,w,2)' will result in the 3x5 array 'out'.

Np=size(A);
if Np(dim)==length(w)
    out=zeros(Np(1:dim-1 dim+1:end));
    for ind1=1:Np(dim)
        out=out+squeeze(A()).*w(ind1);
    end
else disp('Subscripted dimensions dont match.')
end
end
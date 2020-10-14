function list_out = testfcn(list_in)
ind1=1;
% N=size(list_in,1);
list_out=list_in;
while ind1<size(list_out,1)
    indices1=findvec(list_out(1:ind1-1,:),list_out(ind1,:));
    indices2=findvec(list_out(ind1+1:end,:),list_out(ind1,:));
    indices=logical([indices1; 0; indices2]);
    list_out(indices,:)=[];
    ind1=ind1+1;
end
end
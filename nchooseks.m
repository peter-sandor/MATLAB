function out = nchooseks(n,ks)
% This function calculates the multinomial coefficient for a given number of
% tries 'n' and for the # of occurrences for each categories, stored
% as elements of vector 'ks'. The length of vector 'ks' is the number of
% categories (possible outcomes).

Nk=length(ks);
if n==sum(ks)
    denom=1;
    for ind1=1:Nk
        denom=denom*factorial(ks(ind1));
    end
    out=factorial(n)/denom;
else
    disp('n~=sum(ki)!')
    out=[];
end
end
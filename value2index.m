function index = value2index(A,value)
% This function finds the indices of 1D array 'A' that correspond to elements
% of A which have values close to the input 'value'.

Nval = size(value);
temp1 = cumprod(Nval);
Ntot = temp1(end);
% value = reshape(value,[Ntot 1]);
for ind1 = 1:Ntot
    if value(ind1) < min(A) || value(ind1) > max(A)
        index(ind1) = 0;
    else
        temp2 = directfind(A,value(ind1));
        if isempty(temp2)
            smallervals=(A<=value(ind1));
            largervals=(A>=value(ind1));
            if abs(value(ind1)-max(A(logical(smallervals)))) <= abs(min(A(logical(largervals)))-value(ind1))
                index(ind1) = directfind(A,max(A(logical(smallervals))));
            else
                index(ind1) = directfind(A,min(A(logical(largervals))));
            end
        else
            index(ind1) = temp2;
        end
    end
end
if index~=0
    index = reshape(index,Nval);
else
    index = [];
end

    function indx = directfind(B,val)
    % This function finds the elements of 1D array 'B' that have the exact
    % values given by the input 'val', and returns with the indices of
    % these elements.
    smallvals=(B<=val);
    largevals=(B>=val);
    indx=unique(map2colvec(smallvals.*largevals).*map2colvec(1:length(B)));
    indx(indx==0)=[];
    end
    
end
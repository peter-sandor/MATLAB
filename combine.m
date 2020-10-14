function output = combine(A)
% This function takes the cell array 'A' as input, where each cell is a 1D array. Then it creates a 2D
% array 'output' which consists of row vectors containing all possible combinations
% of the elements taken from the input arrays (one element from each).
Np=[];
for ind1=1:length(A)
    Np=[Np length(A{ind1})];
end
NNp=length(Np);
cp=cumprod(Np);
output=zeros([cp(end) NNp]);
for ind1=1:size(output,1)
    temp=[];
    for ind2=1:NNp
        if ind2==1
            cpr=1;
        else cpr=cumprod(Np(1:ind2-1));
        end
        temp=[temp A{ind2}(modP(ceil(ind1/cpr(end)),length(A{ind2})))];
    end
    output(ind1,:)=temp;
end
end
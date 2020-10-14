function out = uniquevec(A)
% This function finds the unique occurrences of each row in array A, and returns with them sorted by the first column.
B=sortrows(A);
out=[];
while ~isempty(B)
    out=[out; B(1,:)];
    indvec=findvec(B,B(1,:));
    B=B(~indvec,:);
end
end
function B = extend(A,N)
% This function creates an array B from A with one extra dimension. The
% size of this extra dimension is N and each element is A along it.
dim=length(size(A));
B=A;
for ind1=1:N-1
	B=cat(dim+1,B,A);
end
B=squeeze(B);
end
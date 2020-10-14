function varargout = findvec(A,x)
% This function finds the vector 'x' amongst the rows of the 2-dim array 'A'.
% number of columns of 'A' should equal to the length of 'x'.
x=map2rowvec(x); % 'x' is transformed to be row vector, if it's not one already.
N2=size(A,2);
if N2==size(x,2)
%     temp=A-ones([size(A,1) 1])*x;
    indvec=ones([size(A,1) 1]);
    for ind1=1:N2
        indvec=indvec.*(A(:,ind1)==x(ind1));
    end
    varargout{1}=logical(indvec);
else varargout{1}=[];
%     'Dimensions dont match.'
end
end
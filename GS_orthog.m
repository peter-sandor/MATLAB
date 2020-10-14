function [Q,R] = GS_orthog(A)
% compute QR using Gram-Schmidt
% [Q, R] = GS_orthog(A) returns a matrix Q with orthonormal columns and an invertible
% upper triangular matrix R so that A = Q*R.
[m,n] = size(A);
for j = 1:n
   v = A(:,j);
   for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
   end
   R(j,j) = norm(v);
   Q(:,j) = v/R(j,j);
end
end
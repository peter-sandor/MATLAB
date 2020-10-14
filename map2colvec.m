function x=map2colvec(x)
% checks if x is row or column vector, and makes a column vector out of it.
% (if input x = [1xN] vector, output will be [Nx1])
[N,M]=size(x);
if M>=2 && N==1
    x=x.';
end
end
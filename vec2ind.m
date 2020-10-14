function xout=vec2ind(x)
% This code returns the indices at which a binary vector x (containing only
% zeros and ones) has a value of 1. If all values are zero, the code
% returns with an empty array ('[]').
N=length(x);
xout=double(map2rowvec(x)).*(1:N); %element by element multiplication
xout(xout==0)=[];
end
function out = sigdig(in,N)
% Finds and returns with the first N significant digits of the input.
alldig=in./10.^floor(log10(abs(in)));
if sign(in)>=0
    out=floor(alldig.*10.^(N-1))./10.^(N-1);
else
    out=ceil(alldig.*10.^(N-1))./10.^(N-1);
end
end


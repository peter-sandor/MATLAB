function out = corrl(in1,in2)

% this code calculates the correlation between two functions:
% h(tau)=integrate f(t)g(t-tau) dt

% f(t)=in1, g(t)=in2, out=h(tau)
% 'in1' and 'in2' have to be the same size, and 'out' will be column vector
% optional argument: supply 1 if 'out' needs to be the full length vector
% / 3*length(in1) /

N2=length(in1);
in1=map2colvec(in1);
in2=[zeros([N2 1]); map2colvec(in2); zeros([N2 1])];
out=zeros([2*N2 1]);

for ind1=1:2*N2
    in2_shifted=circshift(in2,[-N2+ind1 1]);
    out(ind1)=sum(in1.*in2_shifted(N2+1:2*N2));
end

end
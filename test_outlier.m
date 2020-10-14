function [outflag] = test_outlier(array_in)
% Implement Student's t-test to identify a single outlier in a series of
% measurements. 'array_in' is a 1D array.
p=0.999; % confidence level
N0=length(array_in);
non_nans=~isnan(array_in);
index_to_use=vec2ind(non_nans); %gets index of 1 values. 
N1=length(index_to_use);
quantile=tinv(p,N1-1); %What is quantile?
% [tmean,tsigma]=tstatP(N1);
% sigmalevel=3*tsigma;
outflag=zeros([N0 1]);
for ind1=1:N1
    t1(index_to_use(ind1))=student_t(array_in(index_to_use(ind1)),array_in(index_to_use([1:ind1-1 ind1+1:end])));
    if abs(t1(index_to_use(ind1)))<=quantile
        outflag(index_to_use(ind1))=1; % datapoints to keep
    else
        outflag(index_to_use(ind1))=0; % datapoints to discard
    end
end
end
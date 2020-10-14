function corr_out = covar2corr(covar_in)
% calculate the correlation matrix from the covariance matrix
for ind1=1:size(covar_in,1)
    for ind2=1:size(covar_in,2)
        corr_out(ind1,ind2)=covar_in(ind1,ind2)/sqrt(covar_in(ind1,ind1)*covar_in(ind2,ind2));
    end
end
end
function FIT_out = eval_chi2(FIT)

x=FIT.x;
y=FIT.y;
N=length(y);
Nc=length(FIT.param);
for ind1=1:length(FIT.param)
    eval([FIT.param{ind1} '=' num2str(FIT.m(ind1)) ';']);
end    
fit_eq=put_array_ops(FIT.eq);
eval(['y_fit=' fit_eq ';']);

FIT_out=FIT;
FIT_out.chi2=sum((y-y_fit).^2./y);
FIT_out.chi2_reduced=FIT_out.chi2/(N-Nc);
end
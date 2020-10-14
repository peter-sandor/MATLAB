function out = separate_terms(x,y,eqn,paramname,paramval)
% This function takes (x,y) data, the fitting equation and the parameters from the result of a fit
% with the ezfit toolbox, then separates and evaluates the additive terms in the equation. 

N=length(x);
Neqn=length(eqn);
ind_plus = strfind(eqn,'+');
if ~isempty(ind_plus)
    Nterms=length(ind_plus)+1;
    eqn_term{1}=put_array_ops(eqn(1:ind_plus(1)-1));
    for ind1=2:Nterms-1
        eqn_term{ind1}=put_array_ops(eqn(ind_plus(ind1-1)+1:ind_plus(ind1)-1));
    end
    eqn_term{Nterms}=put_array_ops(eqn(ind_plus(Nterms-1)+1:Neqn));
else
    eqn_term{1}=eqn;
    Nterms=1;
end


for ind1=1:Nterms
    eval([paramname{ind1} '=' num2str(paramval(ind1)) ';']);
end
out=zeros([N Nterms]);
for ind1=1:Nterms
    eqn_term{ind1}(strfind(eqn_term{ind1},paramname{ind1}))='1';
    eval(['out(:,' num2str(ind1) ')=' eqn_term{ind1} ';']);
end

end

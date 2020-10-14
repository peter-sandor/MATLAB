clear;
sigma=5;
Nrand=10;
N=100;
x=map2colvec(0:2*pi/(N-1):2*pi);
% y=1 + 1*x + 2*x.^2 + randn([1 N])*sigma;
for ind1=1:N
    y2(ind1,:)=1*ones([1 Nrand]) + 1*x(ind1)*ones([1 Nrand]) + 2*x(ind1).^2*ones([1 Nrand])  + randn([1 Nrand])*sigma;
end
y=mean(y2,2);
sigma2=1/sqrt(Nrand)*map2rowvec(std(y2,[],2));
%%
eqn_terms{1}='A';
eqn_terms{2}='+B*x';
eqn_terms{3}='+C*x^2';
eqn_terms{4}='+D*x^3';
% eqn_terms{5}='+E*x^4';
Nterms=length(eqn_terms);
eqn=[];
for ind1=1:Nterms
    eqn=[eqn eqn_terms{ind1}];
    FIT=ezfit(x,y,['y=' eqn ';A=1;B=1;C=1;D=1;']);
    FIT.dy=sigma2;
    FIT2(ind1)=eval_fit(FIT);
    chi2(ind1)=FIT2(ind1).chi2;
    chi2_red(ind1)=FIT2(ind1).chi2_red;
    if ind1>1
        [Pval(ind1-1),Fval(ind1-1)] = F_test(N,ind1-1,ind1,chi2(ind1-1),chi2(ind1));
    end
end
figure;
subplot(121)
plot(x,y);
title('1 + x + 2*x^2 + noise')
showfit(FIT);
subplot(122)
plot(1:Nterms,chi2_red,'ro')
set(gca,'xtick',1:Nterms)
xlim([2 Nterms])
ylim([0 max(chi2_red(2:Nterms))+1])
xlabel('nr. of terms included')
ylabel('\chi^2_{\nu}')
setfigP
for ind1=2:Nterms
    text(ind1,chi2_red(ind1)+0.5,num2str(Pval(ind1-1)));
end

% str1=[];
% for ind1=1:length(FIT2(Nterms).m)
%     str1=[str1 num2str(sqrt(FIT.B(ind1,ind1))/FIT.m(ind1)) ', '];
% end
% str1=str1(1:end-2);
% disp(str1);
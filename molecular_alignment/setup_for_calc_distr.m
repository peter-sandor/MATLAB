[cos2_M,cos1_M]=calc_cos2(simresult_M.amps0,simresult_M.probs_JM,simresult_M.delay0,simresult_M.delay1,[0.2026 3.46e-8]);
%%
ind_d0=max(vec2ind(simresult_F.delay<0.2));
delay_F0=simresult_F.delay(1:ind_d0)/24.2e-6;
delay_F1=simresult_F.delay(ind_d0+1:end)/24.2e-6;
[cos2_F,cos1_F]=calc_cos2_v2(simresult_F.amps,simresult_M.probs_JM,delay_F0,delay_F1,[0.2026 3.46e-8]);
%%
inF.simresult=simresult_F;
inF.basisnames{1}='sin(K*\theta)^M';
inF.K = [1];
inF.M = [2];
inF=generate_basis(inF);

inM.simresult=simresult_M;
inM.basisnames{1}='cos(K*\theta)^M';
inM.K = [1];
inM.M = [2];
inM=generate_basis(inM);
%%
figure;hold on;
% plot(simresult_M.delay,cos2_M,'k')
% plot(simresult_F.delay,simresult_F.moments(:,1),'b')
% plot((delay_F1-delay_F0(end))*24.2e-6,cos2_F,'r')
plot(simresult_M.delay,inM.basis,'g')
plot(simresult_F.delay,inF.basis,'m')
xlabel('delay [ps]')
ylabel('<cos(\theta)^2>')
legend('MATLAB amps, MATLAB formula','FORTRAN amps, FORTRAN formula','FORTRAN amps, MATLAB formula','MATLAB amps, MATLAB prob. dist.','FORTRAN amps, FORTRAN prob dist.')
%%
calc_prob_vs_theta(ones(size(simresult.amps0)),ones(size(simresult.JM_probs)),simresult.delay0,simresult.delay1(1:10:end),[0.2026 0],simresult.Ylm_red,simresult.theta);
%%
calc_prob_vs_theta(ones(size(out.amps)),ones(size(out.AMP)),out.t_ip,out.TT,[0.2026 0],out.PLeg,map2colvec(out.theta));
%%
calc_prob_vs_theta(simresult0.amps0,simresult0.JM_probs,simresult0.delay0,simresult0.delay1(1:10:end),[0.2026 0],simresult0.Ylm_red,simresult0.theta);
%%
[a,b,c,d]=calc_prob_vs_theta(simresult.amps0,simresult.JM_probs,simresult.delay0,simresult.delay1(1:10:end),[0.2026 0],simresult.Ylm_red,simresult.theta);
%%
[a,b,c,d]=calc_prob_vs_theta(out.amps,out.AMP,out.t_ip,out.TT,[0.2026 0],out.PLeg,map2colvec(out.theta));
%%
figure;
hold on;
plot(a*24.2e-6,c,'k')
plot(a*24.2e-6,d,'r')
ylabel('<cos(\theta)^n>')
xlabel('delay [ps]')
legend('n = 2','n = 4')

figure;imagescP(a*24.2e-6,simresult.theta/pi*180,abs(b).')
xlabel('delay [ps]')
ylabel('\theta [deg]')
title('angular distribution')
%% Check on Spherical Harmonics / Legendre polynomials
Jind=8;
Mind=1;
figure;plot(out.theta,out.PLeg(:,Jind,Mind)/max(out.PLeg(:,Jind,Mind)),'k<')
hold on;plot(simresult.theta,simresult.Ylm_red(:,Jind,Mind)/max(simresult.Ylm_red(:,Jind,Mind)),'ro')
title(['norm ratio = ' num2str(max(out.PLeg(:,Jind,Mind))/max(simresult.Ylm_red(:,Jind,Mind)))])
%%
indJST=16;indM=2;figure;imagescP(abs(squeeze(simresult2.RP(:,:,indM,indJST)+1i*simresult2.IP(:,:,indM,indJST))).^2)
%%
for ind1=1:size(simresult2.RP,3)
    PROB(:,:,ind1)=simresult2.AMP(ind1)*sum(abs(squeeze(simresult2.RP(:,:,1,ind1)+1i*simresult2.IP(:,:,1,ind1))).^2,3) ...
        + 2*simresult2.AMP(ind1)*sum(abs(squeeze(simresult2.RP(:,:,2:end,ind1)+1i*simresult2.IP(:,:,2:end,ind1))).^2,3);
end
figure;imagescP(squeeze(sum(PROB,3)));
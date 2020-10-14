a0=5.29e-11; % [m]
eps0=1/4/pi;
c_au=1/137;
t_scale=24.18884e-18; % [s]
I_scale=3.51e4; % conversion from TW/cm2 to atomic units (1 a.u.=35.1e15 TW/cm2)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
kB=1.38e-23;    % [J/K]
Trot=10;
A=5.097;
% A=0.1;
B=0.4434;
D=8.8e-7;
Jvec=0:30;
Kvec=-30:30;
for indK=1:length(Kvec)
    K=Kvec(indK);
    E_SI(:,indK)=h_SI*c_SI*100*(B*Jvec.*(Jvec+1)+(A-B)*K^2+D*Jvec.^2.*(Jvec+1).^2);
    normfactor0(:,indK) = map2colvec(2.*Jvec+1).*exp(-E_SI(:,indK)/kB/Trot);
end
% figure;imagescP(Kvec,Jvec,normfactor0)
normfactor=sum(sum(normfactor0));
for indK=1:length(Kvec)
    probs_JMK(:,indK)=map2colvec(exp(-E_SI(:,indK)/kB/Trot))/normfactor;
    maxJ(indK)=max(vec2ind(probs_JMK(:,indK)>max(probs_JMK(:,indK)/1000)));
end
% figure;imagescP(Kvec,Jvec,probs_JMK);
max(maxJ)
totProb=sum(sum(probs_JMK.*extend(map2colvec(2.*Jvec+1),2*max(Jvec)+1)));
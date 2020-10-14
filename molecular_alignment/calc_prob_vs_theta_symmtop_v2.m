function varargout = calc_prob_vs_theta_symmtop_v2(amps0,JM_probs,t_ip,t_free,consts,WignerD,theta)
% consts rotational constant 'B' in [1/cm] and centrifugal distortion 'D'
% in [1/cm] with the format: consts = [B D]

fprintf('Calculating Psi(t,theta)')
qE=1.60218e-19;
a0=5.29e-11; % [m]
eps0=1/4/pi;
c_au=1/137;
t_scale=24.188e-18; % [s]
I_scale=3.51e4; % conversion from TW/cm2 to atomic units (1 a.u.=35.1e15 TW/cm2)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
kB=1.38e-23; % [J/K]

Ndelay0=length(t_ip);
Ndelay1=length(t_free);
MaxM=(size(amps0,3)-1);
MaxK=(size(amps0,4)-1);
Mvec=-MaxM:MaxM;
Kvec=-MaxK:MaxK;
MaxJ=MaxM;
maxJ=size(amps0,2)-1;
Jvec=0:maxJ;
amps1=amps0;
amps1(:,:,2:end,2:end,:)=sqrt(2)*amps1(:,:,2:end,2:end,:);
amps1_end=squeeze(amps1(end,:,:,:));
A=consts(2)*100;
B=consts(2)*100;
D=consts(3)*100;
% E_au=map2colvec(h_SI*c_SI/(27.211385*qE)*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2));


for indK=1:length(Kvec)
    K=Kvec(indK);
    E_au(:,indK)=h_SI*c_SI/(27.211385*qE)*(B*Jvec.*(Jvec+1)+(A-B)*K^2+D*Jvec.^2.*(Jvec+1).^2);
%     normfactor = normfactor + sum(map2colvec(2.*Jvec+1).*exp(-E_au(:,indK)/kB/input.Trot));
end
%% moments during interaction with the laser (in case of 'infinite molecules' calculation)
t1=t_ip;
Ntheta=length(theta);
t2=t_free;
Ndelay2=length(t2);
Ndelay=length(t1);
Nmod=100;

% theta2=permute(extend(theta,Ndelay),[2 1]);
D_red=WignerD;
moment1=zeros([Ndelay 1]);
moment2=zeros([Ndelay 1]);
D_sum=zeros([Ndelay Ntheta]);
% Ylm_sum=zeros([Ndelay Ntheta size(amp1,4)]);
D_norm_=zeros([Ndelay Ntheta]);
for ind5=1:Ndelay
    for indM=1:length(Mvec)
        indM2=abs(indM-MaxM-1)+1;
        for indK=1:length(Kvec)
            indK2=abs(indK-MaxK-1)+1;
            for indJSTATE=1:size(amps1,5)
                D_sum(ind5,:)=D_sum(ind5,:) + JM_probs(indJSTATE,indK2)*map2rowvec(sum(squeeze(abs(squeeze(sum(D_red(:,:,indM,indK).*permute(extend(squeeze(amps1(ind5,:,indM2,indK2,indJSTATE)),Ntheta),[2 1]),2))).^2),2));
            end
        end
    end
    D_norm(ind5,:)=D_sum(ind5,:)./sum(D_sum(ind5,:).*map2rowvec(sin(theta)));
	moment1(ind5)=sum(sum(squeeze(D_norm(ind5,:)).*sin(theta).*sin(theta).^2));
    moment2(ind5)=sum(sum(squeeze(D_norm(ind5,:)).*sin(theta).*sin(2*theta).^2));
    if mod(ind5,Nmod)==0
        fprintf('.');
    end
end
%% moments during field-free evolution (in case of 'infinite molecules' approximation)
moment1_ff=zeros([Ndelay2 1]);
moment2_ff=zeros([Ndelay2 1]);
D_sum_ff=zeros([Ndelay2 Ntheta]);
D_norm_ff=zeros([Ndelay2 Ntheta]);
for ind5=1:Ndelay2
    for indK=1:length(Kvec)
        indK2=abs(indK-MaxK-1)+1;
%         E_au_ext=extend(exp(-1i*E_au(:,indK)*(t2(ind5)-t1(end))),size(amps1,3));
        E_au_ext=exp(-1i*E_au(:,indK)*(t2(ind5)-t1(end)));
        for indM=1:length(Mvec)
            indM2=abs(indM-MaxM-1)+1;
            for indJSTATE=1:size(amps1,5)
    %           temp1(ind5,:,:,ind6)=abs(squeeze(sum(Ylm_red.*permute(extend(squeeze(amp1(end,:,:,ind6)).*E_au_ext,Ntheta),[3 1 2]),2))).^2;
                D_sum_ff(ind5,:)=D_sum_ff(ind5,:) + JM_probs(indJSTATE,indK2)*map2rowvec(sum(squeeze(abs(squeeze(sum(D_red(:,:,indM,indK).*permute(extend(squeeze(amps1(end,:,indM2,indK2,indJSTATE)).*E_au_ext,Ntheta),[2 1]),2))).^2),2));
            end
        end
    end
    D_norm_ff(ind5,:)=D_sum_ff(ind5,:)./sum(D_sum_ff(ind5,:).*map2rowvec(sin(theta)));
    moment1_ff(ind5)=sum(sum(squeeze(D_norm_ff(ind5,:)).*sin(theta).*sin(theta).^2));
    moment2_ff(ind5)=sum(sum(squeeze(D_norm_ff(ind5,:)).*sin(theta).*sin(2*theta).^2));
    if mod(ind5,Nmod)==0
        fprintf('.');
    end
end

% varargout{1}=[map2colvec(t1); map2colvec(t2)];
varargout{1}=[D_norm; D_norm_ff];
varargout{2}=[map2colvec(moment1); map2colvec(moment1_ff)];
varargout{3}=[map2colvec(moment2); map2colvec(moment2_ff)];
% varargout{4}=[Ylm_sum; Ylm_sum_ff];
fprintf('done.\n')
end
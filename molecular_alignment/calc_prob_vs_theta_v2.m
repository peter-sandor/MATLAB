function varargout = calc_prob_vs_theta(amps,JM_probs,t_ip,t_free,consts,Ylm,theta)
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
MaxM=size(amps,3)-1;
MaxJ=MaxM;
maxJ=size(amps,2)-1;
Jvec=0:maxJ;
amps_end=squeeze(amps(end,:,:,:));
B=consts(1)*100;
D=consts(2)*100;
E_au=map2colvec(h_SI*c_SI/(27.211385*qE)*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2));

%% moments during interaction with the laser (in case of 'infinite molecules' calculation)
t1=t_ip;
Ntheta=length(theta);
t2=t_free;
Ndelay2=length(t2);
Ndelay=length(t1);
Nmod=100;
amp1=zeros(size(amps));
amp1(:,:,1,:)=amps(:,:,1,:);
amp1(:,:,2:end,:)=sqrt(2)*amps(:,:,2:end,:);
% theta2=permute(extend(theta,Ndelay),[2 1]);
Ylm_red=Ylm;
moment1=zeros([Ndelay 1]);
moment2=zeros([Ndelay 1]);
Ylm_sum=zeros([Ndelay Ntheta]);
% Ylm_sum=zeros([Ndelay Ntheta size(amp1,4)]);
Ylm_norm_=zeros([Ndelay Ntheta]);
for ind5=1:Ndelay
    for ind6=1:size(amp1,4)
        temp0=zeros([Ntheta size(amp1,3)]);
        for ind7=1:Ntheta
            temp0(ind7,:) = squeeze(temp0(ind7,:)) + abs(squeeze(sum(squeeze(Ylm_red(ind7,:,:)).*squeeze(amp1(ind5,:,:,ind6)),1))).^2;
        end
%         temp0=abs(squeeze(sum(temp0,2))).^2;
        Ylm_sum(ind5,:) = Ylm_sum(ind5,:) + JM_probs(ind6)*map2rowvec(sum(squeeze(temp0),2));
    end
    Ylm_norm(ind5,:)=Ylm_sum(ind5,:)./sum(Ylm_sum(ind5,:).*map2rowvec(sin(theta)));
	moment1(ind5)=sum(sum(squeeze(Ylm_norm(ind5,:)).*sin(theta).*sin(theta).^2));
    moment2(ind5)=sum(sum(squeeze(Ylm_norm(ind5,:)).*sin(theta).*sin(2*theta).^2));
    if mod(ind5,Nmod)==0
        fprintf('.');
    end
end
%% moments during field-free evolution (in case of 'infinite molecules' approximation)
moment1_ff=zeros([Ndelay2 1]);
moment2_ff=zeros([Ndelay2 1]);
Ylm_sum_ff=zeros([Ndelay2 Ntheta]);
Ylm_norm_ff=zeros([Ndelay2 Ntheta]);
for ind5=1:Ndelay2
    E_au_ext=extend(exp(-1i*E_au*(t2(ind5)-t2(1))),size(amp1,3));
    for ind6=1:size(amp1,4)
        temp1=zeros([Ntheta size(amp1,3)]);
        for ind7=1:Ntheta
            temp1(ind7,:)=squeeze(temp1(ind7,:)) + abs(squeeze(sum(squeeze(Ylm_red(ind7,:,:)).*squeeze(amp1(end,:,:,ind6)).*E_au_ext,1))).^2;
        end
        if ind5==100
            1;
        end
%         temp1=abs(squeeze(sum(temp1,2))).^2;
%         temp1(ind5,:,:,ind6)=abs(squeeze(sum(Ylm_red.*permute(extend(squeeze(amp1(end,:,:,ind6)).*E_au_ext,Ntheta),[3 1 2]),2))).^2;
        Ylm_sum_ff(ind5,:)=Ylm_sum_ff(ind5,:) + JM_probs(ind6)*map2rowvec(sum(squeeze(temp1),2));
%         disp(num2str(ind6))
    end
    Ylm_norm_ff(ind5,:)=Ylm_sum_ff(ind5,:)./sum(Ylm_sum_ff(ind5,:).*map2rowvec(sin(theta)));
	moment1_ff(ind5)=sum(sum(squeeze(Ylm_norm_ff(ind5,:)).*sin(theta).*sin(theta).^2));
    moment2_ff(ind5)=sum(sum(squeeze(Ylm_norm_ff(ind5,:)).*sin(theta).*sin(2*theta).^2));
    if mod(ind5,Nmod)==0
        fprintf('.');
    end
end

% varargout{1}=[map2colvec(t1); map2colvec(t2)];
varargout{1}=[Ylm_norm; Ylm_norm_ff];
varargout{2}=[map2colvec(moment1); map2colvec(moment1_ff)];
varargout{3}=[map2colvec(moment2); map2colvec(moment2_ff)];
% varargout{4}=[Ylm_sum; Ylm_sum_ff];
fprintf('done.\n')
end
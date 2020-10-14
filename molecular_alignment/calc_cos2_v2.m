function varargout = calc_cos2_v2(amps,probs_JM,t_ip,t_free,consts)
% consts rotational constant 'B' in [1/cm] and centrifugal distortion 'D'
% in [1/cm] with the format: consts = [B D]

fprintf('Calculating <cos^2(theta)>')
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
maxJ=size(amps,1)-1;
Jvec=0:maxJ;
% amps_end=squeeze(amps(end,:,:,:));
amps_end=amps;
[cJJ,cJJp1,cJJp2]= generate_Y_coupling(maxJ);
% [VJJ,VJ1,VJ2]= generate_J_coupling_Bob(maxJ);

B=consts(1)*100;
D=consts(2)*100;
E_au=map2colvec(h_SI*c_SI/(27.211385*qE)*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2));

% cos2_theta0=zeros([Ndelay0 1]);
% cos2_theta1=zeros([Ndelay1 1]);

% for indT=1:Ndelay0
% 	for indM=1:MaxM
%         for indJSTATE=1:MaxJ+1
%             cos2_theta0(indT,indM,indJSTATE)=sum(probs_JM(indJSTATE)*abs(map2colvec(squeeze(amps(indT,:,indM,indJSTATE)))).^2.*(4/3*sqrt(pi/5)*cJJ(:,indM)+1/3)) ...
%             + 8/3*sqrt(pi/5)*probs_JM(indJSTATE)*sum(real(conj(map2colvec(squeeze(amps(indT,1:end-2,indM,indJSTATE)))).*map2colvec(squeeze(amps(indT,3:end,indM,indJSTATE)))).*cJJp2(:,indM)); 
%             cos1_theta0(indT,indM,indJSTATE)= 4*sqrt(pi/3)*probs_JM(indJSTATE)*sum(real(conj(map2colvec(squeeze(amps(indT,1:end-1,indM,indJSTATE)))).*map2colvec(squeeze(amps(indT,2:end,indM,indJSTATE)))).*cJJp1(:,indM)); 
%         
%         end
%         if indM>1
%             cos2_theta0(indT,indM,:)=2*cos2_theta0(indT,indM,:);
%             cos1_theta0(indT,indM,:)=2*cos1_theta0(indT,indM,:);
%         end           
%     end
%     if mod(indT,100)==0
%         fprintf('.')
%     end
% end
% cos2_theta0=squeeze(sum(sum(cos2_theta0,3),2));
% cos1_theta0=squeeze(sum(sum(cos1_theta0,3),2));
cos2_theta0=[];
cos1_theta0=[];

for indTT=1:Ndelay1
	for indM=1:MaxM
        for indJSTATE=1:MaxJ+1
            cos2_theta1(indTT,indM,indJSTATE)=sum(probs_JM(indJSTATE)*abs(squeeze(amps_end(:,indM,indJSTATE))).^2.*(4/3*sqrt(pi/5)*cJJ(:,indM)+1/3)) ...
            + 8/3*sqrt(pi/5)*probs_JM(indJSTATE)*sum(real(conj(squeeze(amps_end(1:end-2,indM,indJSTATE))).*squeeze(amps_end(3:end,indM,indJSTATE)).*map2colvec(exp(-1i*(E_au(3:end)-E_au(1:end-2))*(t_free(indTT)-t_ip(end))))).*cJJp2(:,indM)); 
            cos1_theta1(indTT,indM,indJSTATE)= 4*sqrt(pi/3)*probs_JM(indJSTATE)*sum(real(conj(map2colvec(squeeze(amps_end(1:end-1,indM,indJSTATE)))).*map2colvec(squeeze(amps_end(2:end,indM,indJSTATE))).*map2colvec(exp(-1i*(E_au(2:end)-E_au(1:end-1))*(t_free(indTT)-t_ip(end))))).*cJJp1(:,indM)); 
        end
        if indM>1
            cos2_theta1(indTT,indM,:)=2*cos2_theta1(indTT,indM,:);
            cos1_theta1(indTT,indM,:)=2*cos1_theta1(indTT,indM,:);
        end    
    end
    if mod(indTT+Ndelay0,100)==0
        fprintf('.')
    end
end
cos2_theta1=squeeze(sum(sum(cos2_theta1,3),2));
cos1_theta1=squeeze(sum(sum(cos1_theta1,3),2));

varargout{1}=[map2colvec(cos2_theta0); map2colvec(cos2_theta1)];
varargout{2}=[map2colvec(cos1_theta0); map2colvec(cos1_theta1)];
fprintf('done.\n')
end
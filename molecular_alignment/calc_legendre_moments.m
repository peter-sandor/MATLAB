function varargout = calc_legendre_moments(amps,probs_JM,t_ip,t_free,consts)
% consts rotational constant 'B' in [1/cm] and centrifugal distortion 'D'
% in [1/cm] with the format: consts = [B D]

fprintf('Calculating moments of Legendre polynomials')
qE=1.60218e-19;
a0=5.29e-11; % [m]
eps0=1/4/pi;
c_au=1/137;
t_scale=24.188e-18; % [s]
I_scale=3.51e4; % conversion from TW/cm2 to atomic units (1 a.u.=35.1e15 TW/cm2)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
kB=1.38e-23; % [J/K]

NLgndr=6;
Ndelay0=length(t_ip);
Ndelay1=length(t_free);
MaxM=size(amps,3)-1;
Mvec=0:MaxM;
MaxJ0=MaxM;
MaxJ=size(amps,2)-1;
temp=squeeze(sum(abs(amps(end,:,1,:)),4));
MaxJ1=min(vec2ind(temp<=max(temp)/1000))+1;
if isempty(MaxJ1)
    MaxJ1=MaxJ;
end
Jvec=0:MaxJ1;
kvec=0:2:2*NLgndr-1; % calculate only even order; odd order moments will be zero

B=consts(1)*100;
D=consts(2)*100;
E_au=map2colvec(h_SI*c_SI/(27.211385*qE)*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2));
delta_t = t_free(2)-t_free(1);

moments0a=zeros([Ndelay0 MaxM+1 NLgndr]);
moments1a=zeros([Ndelay1 MaxM+1 NLgndr]);
moments1d=zeros([Ndelay1 MaxM+1 NLgndr NLgndr-1]);
% moments0a=zeros([Ndelay0 MaxJ1+1 MaxJ1+1 MaxM+1 MaxJ+1 NLgndr]);
% moments1a=zeros([Ndelay1 MaxJ1+1 MaxJ1+1 MaxM+1 MaxJ+1 NLgndr]);

for indk=1:NLgndr
    for indM=1:MaxM
        for indJ0=1:MaxJ0+1
            for indJ=1:MaxJ1+1
%                     moments0a(indT,indJ,indJ,indM,indJSTATE,indk) = moments0a(indT,indJ,indJ,indM,indJSTATE,indk) + probs_JM(indJSTATE)*abs(amps(indT,indJ,indM,indJSTATE))^2*(-1)^Mvec(indM)*(2*Jvec(indJ)+1)*sqrt((2*kvec(indk)+1)/4/pi)*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
                moments0a(:,indM,indk) = moments0a(:,indM,indk) + probs_JM(indJ0)*abs(amps(:,indJ,indM,indJ0)).^2*(-1)^Mvec(indM)*(2*Jvec(indJ)+1)*sqrt((2*kvec(indk)+1)/4/pi)*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
%                     moments1a(indTT,indJ,indJ,indM,indJSTATE,indk) = moments1a(indTT,indJ,indJ,indM,indJSTATE,indk) + probs_JM(indJSTATE)*abs(amps(end,indJ,indM,indJSTATE))^2*(-1)^Mvec(indM)*(2*Jvec(indJ)+1)*sqrt((2*kvec(indk)+1)/4/pi)*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
                moments1a(:,indM,indk) = moments1a(:,indM,indk) + probs_JM(indJ0)*abs(amps(end,indJ,indM,indJ0))^2*(-1)^Mvec(indM)*(2*Jvec(indJ)+1)*sqrt((2*kvec(indk)+1)/4/pi)*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJ) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
                for indJp=indJ+1:MaxJ1+1
                    if Jvec(indJp)>=abs(Jvec(indJ)-kvec(indk)) && Jvec(indJp)<=Jvec(indJ)+kvec(indk)
%                       moments0a(indT,indJ,indJp,indM,indJSTATE,indk) = moments0a(indT,indJ,indJp,indM,indJSTATE,indk) + probs_JM(indJSTATE)*2*real(amps(indT,indJ,indM,indJSTATE)*conj(amps(indT,indJp,indM,indJSTATE)))*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]);
                        moments0a(:,indM,indk) = moments0a(:,indM,indk) + probs_JM(indJ0)*2*real(squeeze(amps(:,indJ,indM,indJ0)).*conj(squeeze(amps(:,indJp,indM,indJ0))))*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
%                       moments1a(indTT,indJ,indJ,indM,indJSTATE,indk) = moments1a(indTT,indJ,indJ,indM,indJSTATE,indk) + probs_JM(indJSTATE)*2*real(amps(end,indJ,indM,indJSTATE)*conj(amps(end,indJp,indM,indJSTATE))*exp(-1i*(E_au(indJ)-E_au(indJp))*(t_free(indTT)-t_ip(end))))*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]); 
%                       moments1a(:,indM,indk) = moments1a(:,indM,indk) + probs_JM(indJ0)*2*real(amps(end,indJ,indM,indJ0)*conj(amps(end,indJp,indM,indJ0))*exp(-1i*(E_au(indJ)-E_au(indJp))*(t_free-t_ip(end))))*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]);
                        moments1a(:,indM,indk) = moments1a(:,indM,indk) + probs_JM(indJ0)*2*abs(amps(end,indJ,indM,indJ0)*conj(amps(end,indJp,indM,indJ0)))*cos((E_au(indJp)-E_au(indJ))*(t_free-t_ip(end))+unwrap(angle(amps(end,indJ,indM,indJ0)*conj(amps(end,indJp,indM,indJ0)))))*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]);
                        for indD = 1:NLgndr-1
                            moments1d(:,indM,indk,indD) = moments1d(:,indM,indk,1) + ((E_au(indJp)-E_au(indJ))*delta_t)^indD*probs_JM(indJ0)*2*abs(amps(end,indJ,indM,indJ0)*conj(amps(end,indJp,indM,indJ0)))*cos((E_au(indJp)-E_au(indJ))*(t_free-t_ip(end))+unwrap(angle(amps(end,indJ,indM,indJ0)*conj(amps(end,indJp,indM,indJ0))))+indD*pi/2)*(-1)^Mvec(indM)*sqrt((2*Jvec(indJ)+1)*(2*kvec(indk)+1)*(2*Jvec(indJp)+1)/4/pi)*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[0 0 0])*Wigner3j([Jvec(indJp) kvec(indk) Jvec(indJ)],[-Mvec(indM) 0 Mvec(indM)]);
                        end
                    end
                end
            end
        end
    if indM>1
%             moments0a(indT,:,:,indM,:,indk)=2*moments0a(indT,:,:,indM,:,indk);
        moments0a(:,indM,indk)=2*moments0a(:,indM,indk);
%             moments1a(indTT,:,:,indM,:,indk)=2*moments1a(indTT,:,:,indM,:,indk);
        moments1a(:,indM,indk)=2*moments1a(:,indM,indk);
        moments1d(:,indM,indk,:)=2*moments1d(:,indM,indk,:);
    end
    end
    fprintf('.')
end

% moments0=squeeze(sum(sum(sum(sum(moments0a,5),4),3),2));
moments0=squeeze(sum(moments0a,2));
% moments1=squeeze(sum(sum(sum(sum(moments1a,5),4),3),2));
moments1=squeeze(sum(moments1a,2));
moments1_diff=squeeze(sum(moments1d,2));

varargout{1}=[moments0; moments1];
varargout{2}=[moments1_diff];
fprintf('done.\n')
end
function varargout = solve_align_symmtop(varargin)
% use atomic units, except for some of the constants below
tic;

a0=5.29e-11; % [m]
eps0=1/4/pi;
c_au=1/137;
t_scale=24.18884e-18; % [s]
I_scale=3.51e4; % conversion from TW/cm2 to atomic units (1 a.u.=35.1e15 TW/cm2)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
kB=1.38e-23; % [J/K]

input=varargin{1};
% conversion to atomic units
input.laser_int1=varargin{1}.laser_int1;%/I_scale; % laser peak intensity in [atomic units]
input.laser_int2=varargin{1}.laser_int2;%/I_scale; % laser peak intensity in [atomic units]
input.laser_fwhm1=varargin{1}.laser_fwhm1/t_scale/1e12; % laser pulse duration intensity FWHM in [atomic u.]
input.laser_fwhm2=varargin{1}.laser_fwhm2/t_scale/1e12; % laser pulse duration intensity FWHM in [atomic u.]
input.rot_const_A=varargin{1}.rot_cnst_A*a0*100; % rotational constant 'A' in [a0]
input.rot_const_B=varargin{1}.rot_cnst_B*a0*100; % rotational constant 'B' in [a0]
input.centrif=varargin{1}.centrif*a0*100; % centrifugal distortion 'D' [a0]
input.t0=varargin{1}.t0/t_scale/1e12; % stepsize in [atomic u.];
input.maxdelay=varargin{1}.maxdelay/t_scale/1e12; % max timedelay in [atomic u.]
input.timestep=varargin{1}.timestep/t_scale/1e12; % stepsize in [atomic u.];
sigma1=2*input.laser_fwhm1/(2*sqrt(2*log(2)));
sigma2=2*input.laser_fwhm2/(2*sqrt(2*log(2)));

Jvec=0:input.maxJ;
[cJJ,cJJp1,cJJp2]= generate_D_coupling(input.maxJ);
cJJ=cJJ*pi/3*sqrt(2/5);
cJJp1=cJJp1*2*pi*sqrt(2/3);
cJJp2=cJJp2*pi/3*sqrt(2/5);

Tlimit=2*input.laser_fwhm2;
Ndelay0=100;
t_ip=-Tlimit:(2*Tlimit)/(Ndelay0-1):Tlimit;
t_free=map2colvec(Tlimit+input.timestep:input.timestep:input.maxdelay);
Ndelay1=length(t_free);
% F1=sqrt(2/8.85e-12/3e8*input.laser_int1*1e16)/5.14e11;
% F2=sqrt(2/8.85e-12/3e8*in.laser_int2*1e16)/5.14e11;
F1=sqrt(input.laser_int1/I_scale);
F2=sqrt(input.laser_int2/I_scale);
Ntheta=input.Ntheta;
theta=0:pi/(Ntheta-1):pi;

% for this part only, use SI units
A=input.rot_const_A/a0/100; % [1/cm]
B=input.rot_const_B/a0/100; % [1/cm]
D=input.centrif/a0/100; % [1/cm]
E_SI0=h_SI*c_SI*100*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2);
probs0_JMK=exp(-E_SI0/kB/input.Trot)/sum(exp(-E_SI0/kB/input.Trot));
MaxJ=max(vec2ind(probs0_JMK>max(probs0_JMK)/1000))+1;
MaxM=MaxJ;
MaxK=MaxJ;
Mvec=-MaxM:MaxM;
Kvec=Mvec;
cJJ=cJJ(:,input.maxJ-MaxM+1:input.maxJ+MaxM+1,input.maxJ-MaxK+1:input.maxJ+MaxK+1);
cJJp1=cJJp1(:,input.maxJ-MaxM+1:input.maxJ+MaxM+1,input.maxJ-MaxK+1:input.maxJ+MaxK+1);
cJJp2=cJJp2(:,input.maxJ-MaxM+1:input.maxJ+MaxM+1,input.maxJ-MaxK+1:input.maxJ+MaxK+1);
normfactor=0;
for indK=1:length(Kvec)
    K=Kvec(indK);
    E_SI(:,indK)=h_SI*c_SI*100*(B*Jvec.*(Jvec+1)+(A-B)*K^2+D*Jvec.^2.*(Jvec+1).^2);
    normfactor = normfactor + sum(map2colvec(2.*Jvec+1).*exp(-E_SI(:,indK)/kB/input.Trot));
end
for indK=1:length(Kvec)
    probs_JMK(:,indK)=map2colvec(exp(-E_SI(:,indK)/kB/input.Trot))/normfactor;
end
probs_JK=probs_JMK.*extend((0:input.maxJ)*2+1,2*MaxJ+1);
% end of SI part

coeffs=probs_JK;
ang_distr=zeros([Ntheta Ndelay0]);
ang_distr_free=zeros([Ntheta Ndelay1]);
Jamps=zeros([input.maxJ+1 1]);
cJmk0=zeros([input.maxJ+1, 2*MaxM+1, 2*MaxK+1]);
cJm_end=zeros([input.maxJ+1, 2*MaxM+1, 2*MaxK+1]);
if input.solvetype==1
    cJm_ip=zeros([Ndelay0, input.maxJ+1, 2*MaxM+1, 2*MaxK+1, MaxJ+1]);
    cJm_amp=cJm_ip;
else
    cJm_ip=zeros([Ndelay0, input.maxJ+1, 2*MaxM+1]);
    cJm_amp=cJm_ip;
end
cos_sq_theta=zeros([Ndelay1 1]);
cJm_free{1}=zeros([length(t_free) input.maxJ+1]);
E_au=map2colvec(E_SI/(27.211385*1.60218e-19));

PREC_ABS = 5E-6; %Absolute Precision
PREC_REL = 5E-8; %Relative Precision
options = odeset('RelTol',1E-8,'AbsTol',1E-10,'InitialStep',25,'MaxStep',50); %Options for integrator
fprintf('Solving TDSE')
for indM=1:2*MaxM+1
    Mval=Mvec(indM);
    if Mval==0
        1;
    end
	for indK=1:2*MaxK+1
       Kval=Mvec(indK);
        if input.solvetype==0
    %         cJm0(:,indM)=sqrt(coeffsJ./(2*Jvec+1)).*exp(1i*input.rand*rand([1 input.maxJ+1])*2*pi);
    %         if indM==in.maxJ+1 % initialize only two |J,M> states
    %             indJ0a=8;
    %             indJ0b=10;
    %             cJm0(:,indM)=[zeros([indJ0a 1]); 0.5; zeros([indJ0b-indJ0a-1 1]); 0.5; zeros([in.maxJ-indJ0b 1])];
    %         end
    %         cJm0(abs(Mvec(indM))>Jvec,indM)=0;
    %         [t_out{indM},cJm_out{indM}] = ode45(@Hamiltonian,[-Tlimit Tlimit],cJm0(:,indM),options);
    %         Jamps=Jamps+map2colvec(cJm_out{indM}(end,:));
    %         cJm_end(:,indM)=cJm_out{indM}(end,:);
    %         cJm_ip(:,:,indM)=interp1(t_out{indM},cJm_out{indM},t_ip);
    %         if indM>1
    %             cJm_amp=2*squeeze(cjM_ip(:,:,indM)).*permute(extend(sqrt(coeffsJ./(2*Jvec+1)),Ndelay0),[2 1]);
    %         else
    %             cJm_amp=squeeze(cjM_ip(:,:,indM)).*permute(extend(sqrt(coeffsJ./(2*Jvec+1)),Ndelay0),[2 1]);
    %         end

        elseif input.solvetype==1
            for indJ=1:MaxJ
                Jval=Jvec(indJ);
                cJ_initial=zeros([input.maxJ+1 1]);
                if Jvec(indJ)>=abs(Mval) && Jvec(indJ)>=abs(Kval)
                    cJ_initial(indJ)=exp(1i*input.rand*rand*2*pi);
                    cJmk0(indJ,indM,indK)=cJmk0(indJ,indM,indK) + cJ_initial(indJ);
                    cJmk0(abs(Mval)>Jvec,indM,indK)=0;
                    cJmk0(abs(Kval)>Jvec,indM,indK)=0;
                    [t_out{indJ,indM,indK},cJmk_out{indJ,indM,indK}] = ode45(@Hamiltonian,[-Tlimit Tlimit],cJ_initial,options);
                    Jamps=Jamps+map2colvec(cJmk_out{indJ,indM,indK}(end,:));
                    cJm_ip(:,:,indM,indK,indJ)=interp1(t_out{indJ,indM,indK},cJmk_out{indJ,indM,indK},t_ip);
                end
            end
        end
    %     if mod(indM,5)==0 || indM==MaxM
    %         disp([num2str(indM) '/' num2str(MaxM)]);
    %     end
    end
fprintf('.')
end
cJm_end=squeeze(cJm_ip(end,:,:,:,:));
fprintf('done.\n')
if input.calc_cos2==1
%     [cos2_theta,cos1_theta]=calc_cos2(cJm_ip,probs_JMK,t_ip,t_free,[varargin{1}.rot_cnst_B varargin{1}.centrif]);
else
    cos1_theta=[];
    cos2_theta=[];
end
if input.calc_prob==1
    WignerD=wigner_d_P(theta,input.maxJ);
    WignerD=WignerD(:,:,input.maxJ+1-MaxJ:input.maxJ+1+MaxJ,input.maxJ+1-MaxJ:input.maxJ+1+MaxJ);
    [ang_prob,moment1,moment2]=calc_prob_vs_theta_symmtop(cJm_ip,probs_JMK,t_ip,t_free,[varargin{1}.rot_cnst_A varargin{1}.rot_cnst_B varargin{1}.centrif],WignerD,theta);
else
    WignerD=[];
    ang_prob=[];
    moment1=[];
    moment2=[];
end
out.runtime=toc;
out.input=input;
out.theta=theta;
out.WignerD=WignerD;
out.delay0=t_ip;
out.delay1=t_free;
out.delay_au=[map2colvec(t_ip); map2colvec(t_free)];
out.delay=out.delay_au*t_scale*1e12;
out.probs_JK=probs_JK;
out.probs_JMK=probs_JMK;
out.laser_int=(F2*exp(-t_ip.^2/sigma2^2)).^2;
% out.amps=cJm_amp;
out.amps0=cJm_ip;
out.cos2_theta=cos2_theta;
out.cos1_theta=cos1_theta;
out.moment1=moment1;
out.moment2=moment2;
out.ang_prob=ang_prob;
varargout{1}=out;
       
function cJmk_dot = Hamiltonian(t,cJmk) % include orientation field
    field1=F1*exp(-(t-input.t0).^2/sigma1^2);
    field2=F2*exp(-t.^2/sigma2^2);
    kin_nrg=map2colvec(E_au);
    cJmk_dot=-1i*(squeeze(kin_nrg(:,indK)) - input.polar_antr/12*field2.^2 - input.polar_antr*field2.^2*cJJ(:,indM,indK)).*cJmk ... % Jp = J
        + 1i*input.dipole*field1*[cJJp1(:,indM,indK); 0].*[cJmk(2:end); 0] ... % Jp = J+1
        + 1i*input.polar_antr*field2.^2*[cJJp2(:,indM,indK); 0; 0].*[cJmk(3:end); 0; 0]; % Jp = J+2
    cJmk_dot(2:end) = cJmk_dot(2:end) + 1i*2*pi*sqrt(2/3)*input.dipole*field1*cJJp1(:,indM,indK).*cJmk(1:end-1); % Jp = J-1
    cJmk_dot(3:end) = cJmk_dot(3:end) + 1i*input.polar_antr*field2.^2*cJJp2(:,indM,indK).*cJmk(1:end-2); % Jp = J-2
end

end
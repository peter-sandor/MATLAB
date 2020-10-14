function varargout = solve_align_linear(varargin)
% This code solves the TDSE to calculate molecular alignment and THz orientation, i.e.
% interaction of a short laser pulse + THz pulse with the rotational wavefunction of a
% linear molecule. Calculations with the short pulses are tested and do
% work. Calculations involving the THz field are not tested!
tic;
% use atomic units, except for some of the constants below
a0=5.29e-11; % [m]
eps0=1/4/pi;
c_au=1/137;
t_scale=24.18884e-18; % [s]
I_scale=3.51e4; % conversion from TW/cm2 to atomic units (1 a.u.=35.1e15 TW/cm2)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
kB=1.38e-23; % [J/K]
qE=1.6022e-19; % [C]
E_Ryd=27.211385; % [eV]

input=varargin{1};
% conversion to atomic units
input.laser_int1=varargin{1}.laser_int1;%/I_scale; % laser peak intensity in [atomic units]
input.laser_int2=varargin{1}.laser_int2;%/I_scale; % laser peak intensity in [atomic units]
input.laser_fwhm1=varargin{1}.laser_fwhm1/t_scale/1e12; % laser pulse duration intensity FWHM in [atomic u.]
input.laser_fwhm2=varargin{1}.laser_fwhm2/t_scale/1e12; % laser pulse duration intensity FWHM in [atomic u.]
input.rot_const_B=varargin{1}.rot_const_B*a0*100; % rotational constant 'B' in [a0]
input.centrif=varargin{1}.centrif*a0*100; % centrifugal distortion 'D' [a0]
input.t0=varargin{1}.t0/t_scale/1e12; % stepsize in [atomic u.];
input.maxdelay=varargin{1}.maxdelay/t_scale/1e12; % max timedelay in [atomic u.]
input.timestep=varargin{1}.timestep/t_scale/1e12; % stepsize in [atomic u.];
sigma1=2*input.laser_fwhm1/(2*sqrt(2*log(2)));
sigma2=2*input.laser_fwhm2/(2*sqrt(2*log(2)));

Jvec=0:input.maxJ;
% [cJJ,cJJp1,cJJp2]= generate_J_coupling_Bob(in.maxJ); cJJp2(end-1:end,:)=[]; cJJp2(:,[1 end])=[]; cJJ=abs(cJJ); cJJp2=abs(cJJp2);
[cJJ,cJJp1,cJJp2] = generate_Y_coupling(input.maxJ); % calculate array with coupling coefficients
cJJ=cJJ/3*sqrt(pi/5);
cJJp1=cJJp1*2*sqrt(pi/3);
cJJp2=cJJp2/3*sqrt(pi/5);

Tlimit=2*input.laser_fwhm2;
Ndelay0=100;
t_ip=-Tlimit:(2*Tlimit)/(Ndelay0-1):Tlimit; % generate delay array for times when interaction with laser pulse is not negligible
t_free=map2colvec(Tlimit+input.timestep:input.timestep:input.maxdelay); % generate delay array for times of field-free evolution
Ndelay1=length(t_free);
F1=sqrt(input.laser_int1/I_scale); % THz field amplitude [atomic units]
F2=sqrt(input.laser_int2/I_scale); % laser pulse field amplitude [atomic units]
Ntheta=input.Ntheta;
theta=0:pi/(Ntheta-1):pi; % array for angle theta

% for this part only, use SI units
B=input.rot_const_B/a0/100; % [1/cm]
D=input.centrif/a0/100; % [1/cm]
E_SI=h_SI*c_SI*100*(B*Jvec.*(Jvec+1)+D*Jvec.^2.*(Jvec+1).^2);
E_au=map2colvec(E_SI/(E_Ryd*qE));
% generate thermal distribution of J-states (with 2*J+1 multiplicity)
probs_J=map2rowvec(2.*Jvec+1).*map2rowvec(exp(-E_SI/kB/input.Trot))/sum((2.*Jvec+1).*exp(-E_SI/kB/input.Trot));
% end of SI part

MaxJ=max(vec2ind(probs_J>max(probs_J)/1000))+1;
MaxM=MaxJ;
Mvec=0:MaxM;
% probs_J=[map2rowvec(probs_J5) 0];
probs_J(MaxJ+1:end)=0;
probs_J=probs_J/sum(probs_J); % probabilities with all M-states added up for a given J
probs_JM=probs_J./(Jvec*2+1); % probabilities for a given |J,M> state

% generate reduced spherical harmonics (only theta-dependence, not phi)
legendre_poly=LegendreP(cos(theta),input.maxJ);
Ylm_reduced=zeros([Ntheta input.maxJ+1 MaxM+1]);
for J=0:input.maxJ
    maxM=min([MaxM,J]);
    for M=0:maxM
        Ylm_reduced(:,J+1,M+1)=sqrt((2*J+1)*factorial(J-M)/4/pi/factorial(J+M))*legendre_poly(:,J+1,input.maxJ+1+M);
    end
end

% initialize arrays for storing state amplitudes
ang_distr=zeros([Ntheta Ndelay0]);
ang_distr_free=zeros([Ntheta Ndelay1]);
Jamps=zeros([input.maxJ+1 1]);
cJm0=zeros([input.maxJ+1, MaxM+1]);
cJm_end=zeros([input.maxJ+1, MaxM+1]);
if input.solvetype==1
    cJm_ip=zeros([Ndelay0, input.maxJ+1, MaxM+1, MaxJ+1]);
    cJm_amp=cJm_ip;
else
    cJm_ip=zeros([Ndelay0, input.maxJ+1, MaxM+1]);
    cJm_amp=cJm_ip;
end
cos_sq_theta=zeros([Ndelay1 1]);
cJm_free{1}=zeros([length(t_free) input.maxJ+1]);

options = odeset('RelTol',1E-10,'AbsTol',1E-12,'InitialStep',25,'MaxStep',50); %Options for integrator
fprintf('Solving TDSE')
for indM=1:MaxM+1 % loop for M quantum numbers
    if input.solvetype==0 % this case is obsolete
        % initialize state amplitudes
        cJm0(:,indM)=sqrt(probs_J./(2*Jvec+1)).*exp(1i*input.rand*rand([1 input.maxJ+1])*2*pi);
%         if indM==in.maxJ+1 % initialize only two |J,M> states
%             indJ0a=8;
%             indJ0b=10;
%             cJm0(:,indM)=[zeros([indJ0a 1]); 0.5; zeros([indJ0b-indJ0a-1 1]); 0.5; zeros([in.maxJ-indJ0b 1])];
%         end
        cJm0(abs(Mvec(indM))>Jvec,indM)=0;
        [t_out{indM},cJm_out{indM}] = ode45(@Hamiltonian,[-Tlimit Tlimit],cJm0(:,indM),options);
        Jamps=Jamps+map2colvec(cJm_out{indM}(end,:));
        cJm_end(:,indM)=cJm_out{indM}(end,:);
        cJm_ip(:,:,indM)=interp1(t_out{indM},cJm_out{indM},t_ip); % interpolate
        if indM>1
            cJm_amp=2*squeeze(cjM_ip(:,:,indM)).*permute(extend(sqrt(probs_J./(2*Jvec+1)),Ndelay0),[2 1]);
        else
            cJm_amp=squeeze(cjM_ip(:,:,indM)).*permute(extend(sqrt(probs_J./(2*Jvec+1)),Ndelay0),[2 1]);
        end
        
    elseif input.solvetype==1
        for indJ0=1:MaxJ
            Jval=Jvec(indJ0);
            Mval=Mvec(indM);
            cJ_initial=zeros([input.maxJ+1 1]);
            if Jvec(indJ0)>=abs(Mvec(indM)) % loop for initialized J0 states
                cJ_initial(indJ0)=exp(1i*input.rand*rand*2*pi); % initialize state amplitudes
                cJm0(indJ0,indM)=cJm0(indJ0,indM) + cJ_initial(indJ0);
                cJm0(abs(Mvec(indM))>Jvec,indM)=0;
                % solution is calculated here for a given M and J0, usding ode45
                [t_out{indJ0,indM},cJm_out{indJ0,indM}] = ode45(@Hamiltonian,[-Tlimit Tlimit],cJ_initial,options);
                Jamps=Jamps+map2colvec(cJm_out{indJ0,indM}(end,:));
                % interpolate obtained amplitudes 'cJm' for times where interaction with
                % laser pulse is significant
                cJm_ip(:,:,indM,indJ0)=interp1(t_out{indJ0,indM},cJm_out{indJ0,indM},t_ip);
%                 cJm_ip(:,6:end,indM,indJ0)=0; % this line is for
%                 troubleshooting purposes, normally leave it commented
%                 out
                maxsteptaken(indJ0,indM)=max(diff(t_out{indJ0,indM}));
                minsteptaken(indJ0,indM)=min(diff(t_out{indJ0,indM}));
                meansteptaken(indJ0,indM)=mean(diff(t_out{indJ0,indM}));
            end
        end
    end
%     if mod(indM,5)==0 || indM==MaxM
%         disp([num2str(indM) '/' num2str(MaxM)]);
%     end
fprintf('.')
end
fprintf('done.\n')
runtime1=toc;
tic;
if input.calc_moments==1
    % calculate cos(theta) and cos^2(theta) moments based on analytical formula
%     [cos2_theta,cos1_theta]=calc_cos2(cJm_ip,probs_JM,t_ip,t_free,[varargin{1}.rot_const_B varargin{1}.centrif]);
    [moments moments_diff]=calc_legendre_moments(cJm_ip,probs_JM,t_ip,t_free,[varargin{1}.rot_const_B varargin{1}.centrif]);
else
    cos1_theta=[];
    cos2_theta=[];
    moments=0;
end
runtime2=toc;
tic;
if input.calc_prob==1
    % calculate axis distribution
    [ang_prob,moment1,moment2]=calc_prob_vs_theta_linear(cJm_ip,probs_JM,t_ip,t_free,[varargin{1}.rot_const_B varargin{1}.centrif],Ylm_reduced,theta);
else
    ang_prob=[];
    moment1=[];
    moment2=[];
end
runtime3=toc;
out.runtime=[runtime1 runtime2 runtime3];
out.maxstep=maxsteptaken;
out.minstep=minsteptaken;
out.meanstep=meansteptaken;
out.input=input;
out.theta=theta;
out.delay0=t_ip;
out.delay1=t_free;
out.delay_au=[map2colvec(t_ip); map2colvec(t_free)];
out.delay=out.delay_au*t_scale*1e12;
out.probs_J=probs_J;
out.probs_JM=probs_JM;
out.laser_int=(F2*exp(-t_ip.^2/sigma2^2)).^2;
% out.amps=cJm_amp;
out.amps0=cJm_ip;
% out.cos2_theta=cos2_theta;
% out.cos1_theta=cos1_theta;
out.moment1=moment1;
out.moment2=moment2;
out.moments=moments;
out.moments_diff=moments_diff;
out.ang_prob=ang_prob;
out.Ylm_red=Ylm_reduced;
varargout{1}=out;

% This is the function the ode45 solver uses
function cJmdot = Hamiltonian(t,cJm) % include orientation field
    field1=F1*exp(-(t-input.t0).^2/sigma1^2);
    field2=F2*exp(-t.^2/sigma2^2);
    kin_nrg=map2colvec(E_au);
    cJmdot=-1i*(kin_nrg - input.polar_antr/12*field2.^2 - input.polar_antr*field2.^2*cJJ(:,indM)).*cJm ... % Jp = J
        + 1i*input.dipole*field1*[cJJp1(:,indM); 0].*[cJm(2:end); 0] ... % Jp = J+1
        + 1i*input.polar_antr*field2.^2*[cJJp2(:,indM); 0; 0].*[cJm(3:end); 0; 0]; % Jp = J+2
    cJmdot(2:end) = cJmdot(2:end) + 1i*input.dipole*field1*cJJp1(:,indM).*cJm(1:end-1); % Jp = J-1
    cJmdot(3:end) = cJmdot(3:end) + 1i*input.polar_antr*field2.^2*cJJp2(:,indM).*cJm(1:end-2); % Jp = J-2
end

end
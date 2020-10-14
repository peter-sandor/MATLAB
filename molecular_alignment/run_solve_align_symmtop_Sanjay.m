a0=5.29e-11; % in [m]
eps0=1/4/pi;
c=1/137;
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
qE=1.60218e-19; % [C]
E_Hartree=27.211385; % [eV]

in.laser_int1=0; % laser peak intensity in [TW/cm2]
in.laser_int2=50; % laser peak intensity in [TW/cm2]
in.Trot=5; % rotational temperature in [K]
in.t0=0; % relative delay of orientation pulse with respect to alignment pulse
in.laser_fwhm1=1; % laser pulse duration intensity FWHM in [ps]
in.laser_fwhm2=0.075; % laser pulse duration intensity FWHM in [ps]
in.polar_antr=11.5; % polarizability anisotropy [atomic units]
in.dipole=0.748; % molecular dipole moment [atomic units]
in.rot_cnst_A=1*5.097; % rotational constant 'A' in [1/cm]
in.rot_cnst_B=0.4434; % rotational constant 'B' in [1/cm]
in.centrif=8.8e-7; % centrifugal distortion 'D' [1/cm]
in.abund_evenJ=1; % abundance of even J states
in.abund_oddJ=1; % abundance of odd J states
in.maxJ=40; % max J level considered
in.maxdelay=50; % max timedelay in [ps]
in.timestep=0.1;
in.Ntheta=200;
in.rand=0; % random phase for different (J,M) states? (1: yes, 0: no)
in.solvetype=1; % 0: populate all initial states and solve ('single molecule' calculation); 1: populate a single initial state, solve, then repeat ('infinite number of molecules' approximation)
in.calc_cos2=0; % calculate cos^2(theta)? (1=yes,0=no)
in.calc_prob=0; % calculate angular probability distribution? (1=yes,0=no)

Nrep=1;
Nfree=in.maxdelay/in.timestep;
cos_sq=zeros([Nfree,Nrep]);
% amps=zeros([50 in.maxJ+1 2*in.maxJ+1]);
clear simresult
for ind1=1:Nrep
    simresult(ind1) = solve_align_symmtop_Sanjay(in);
    rotational_distribution(simresult(ind1));
end
%%
figure;
hold on;
plot(simresult.delay,simresult.moment1,'k')
plot(simresult.delay,simresult.moment2,'r')
% ylabel('<cos(\theta)^n>')
xlabel('delay [ps]')
legend('<sin^2(\theta)>','<sin^2(2*\theta)>')

figure;imagescP(simresult.delay,simresult.theta/pi*180,abs(simresult.ang_prob).')
xlabel('delay [ps]')
ylabel('\theta [deg]')
title('angular distribution')
%%
% simresult=simresult_array(11,9);
MaxM=size(simresult.amps0,3);
MaxK=size(simresult.amps0,4);
Jvec=map2colvec(0:simresult.input.maxJ);
E_SI=h_SI*c_SI*100*(simresult.input.rot_const_B*Jvec.*(Jvec+1)+simresult.input.centrif*Jvec.^2.*(Jvec+1).^2);
E_au=map2rowvec(E_SI/(E_Hartree*qE));

Jpop0=zeros([simresult.input.maxJ+1 1]);
Jpop1=zeros([simresult.input.maxJ+1 1]);
Jpop=zeros([length(simresult.delay0) simresult.input.maxJ+1 1]);
amps0=simresult.amps0;
clear amps1 ReAmp phase1;
if length(size(amps0))==5
    for indK=1:MaxK
        for indM=1:MaxM;
            amps1(:,:,indM,indK,:)=squeeze(amps0(:,:,indM,indK,:)).*permute(extend(extend(sqrt(squeeze(simresult.probs_JMK(1:size(amps0,5),indK))),length(simresult.delay0)),size(amps0,2)),[2 3 1]);
            ReAmp1(indM,indK,:)=real(squeeze(sum(amps1(end,:,indM,indK,:),5)));
            phase1(indM,indK,:)=unwrap(angle(squeeze(sum(amps1(end,:,indM,indK,:),5))));
            Jpop1=Jpop1+map2colvec(abs(squeeze(sum(amps1(end,:,indM,indK,:),5))).^2);
            Jpop0=Jpop0+map2colvec(abs(sum(amps1(1,:,indM,indK,:),5)).^2);
            Jpop=Jpop+map2colvec(abs(sum(amps1(:,:,indM,indK,:),5)).^2);
        end
    end
elseif length(size(amps1))==4
    for indM=1:MaxM;
        ReAmp1(indM,:)=real(squeeze(sum(amps(end,:,indM),3)));
        phase1(indM,:)=unwrap(angle(squeeze(sum(amps(end,:,indM),3))));
        Jpop1=Jpop1+map2colvec(abs(squeeze(sum(amps(end,:,indM),3))).^2);
        Jpop0=Jpop1+map2colvec(abs(sum(amps(1,:,indM),3)).^2);
    end
end
% Jpop=Jpop/sum(Jpop);
figure;
% subplot(121)
hold on;
bar(Jvec,Jpop0,'k')
bar(Jvec,Jpop1,'r')
xlabel('J')
legend('before pulse','after pulse');
ylabel('state population');
xlim([0 simresult.input.maxJ])

figure;
plot(simresult.delay0*24e-6,sum(Jpop,2),'k')
xlabel('delay [ps]')
ylabel('total population')
setfigP
a0=5.29e-11; % in [m]
eps0=1/4/pi;
c=1/137;
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
qE=1.6022e-19; % [C]

in.laser_int1=0*1e-1; % laser peak intensity in [TW/cm2]
in.laser_int2=5; % laser peak intensity in [TW/cm2]
in.Trot=5; % rotational temperature in [K]
in.t0=0; % relative delay of orientation pulse with respect to alignment pulse
in.laser_fwhm1=1; % laser pulse duration intensity FWHM in [ps]
in.laser_fwhm2=0.075; % laser pulse duration intensity FWHM in [ps]
in.polar_antr=31; % polarizability anisotropy [atomic units]
in.dipole=0.33; % molecular dipole moment [atomic units]
in.rot_const_B=0.2026; % rotational constant 'B' in [1/cm]
in.centrif=3.46e-8; % centrifugal distortion 'D' [1/cm]
in.abund_evenJ=1; % abundance of even J states
in.abund_oddJ=1; % abundance of odd J states
in.maxJ=40; % max J level considered
in.maxdelay=90; % max timedelay in [ps]
in.timestep=0.1;
in.Ntheta=200;
in.rand=0; % random phase for different (J,M) states? (1: yes, 0: no)
in.solvetype=1; % 0: populate all initial states and solve ('single molecule' calculation); 1: populate a single initial state, solve, then repeat ('infinite number of molecules' approximation)
in.calc_moments=1; % calculate cos^2(theta)? (1=yes,0=no)
in.calc_prob=1; % calculate angular probability distribution? (1=yes,0=no)

Nrep=1;
Nfree=in.maxdelay/in.timestep;
cos_sq=zeros([Nfree,Nrep]);
% amps=zeros([50 in.maxJ+1 2*in.maxJ+1]);
clear simresult;
for ind1=1:Nrep
    simresult(ind1) = solve_align_linear(in);
    rotational_distribution(simresult(ind1));
end
%%
NLgndr=size(simresult.moments,2);
Nbasis=floor(NLgndr);
in2.simresult=simresult;
in2.basisnames{1}='1';
if Nbasis>1
    for ind1=2:Nbasis
        in2.basisnames{ind1}='sin(K*\theta)^M';
    end
end
in2.K = [0:Nbasis-1];
in2.M = [0 2*ones([1 Nbasis-1])];
in2 = generate_basis(in2);

Jvec=0:in.maxJ;
theta=map2colvec(0:pi/(in.Ntheta-1):pi); % array for angle theta
legendre_poly=LegendreP(cos(theta),in.maxJ);
Yl0=zeros([in.Ntheta in.maxJ+1]);
for indJ=1:in.maxJ+1
    M=0;
    J=Jvec(indJ);
    Yl0(:,indJ)=sqrt((2*J+1)*factorial(J-M)/4/pi/factorial(J+M))*legendre_poly(:,indJ,in.maxJ+1+M);
end
norm1=2./(2*Jvec+1);%./((2*Jvec+1)/4/pi);
basis_nonorth_theta = zeros([in.Ntheta Nbasis]);
basis_nonorth_delay = zeros([length(simresult.delay) Nbasis]);
coeff = zeros([NLgndr Nbasis]);
for ind1=1:Nbasis
    if ind1>1
        basis_nonorth_theta(:,ind1)=sin((ind1-1)*theta).^2;
    else
        basis_nonorth_theta(:,1)=ones([in.Ntheta 1]);
    end
%     str_legend{ind1}=['sin^2(' num2str((ind1-1)) '\theta)'];
    for indJ=1:NLgndr
        coeff(indJ,ind1)=2*pi*trapz(theta,basis_nonorth_theta(:,ind1).*Yl0(:,2*indJ-1).*sin(theta));
        basis_nonorth_delay(:,ind1) = basis_nonorth_delay(:,ind1) + coeff(indJ,ind1)*simresult.moments(:,indJ);
    end
end

figure;
hold on;
plot(simresult.delay,in2.basis,'-')
plot(simresult.delay,basis_nonorth_delay,'--')
xlabel('delay [ps]')
str_legend{1}='1';
for ind1=2:Nbasis
    str_legend{ind1}=['sin^2(' num2str((ind1-1)) '\theta)'];
end
legend(cat(2,str_legend,str_legend))
%%
figure;imagescP(simresult.delay,simresult.theta/pi*180,abs(simresult.ang_prob).')
xlabel('delay [ps]')
ylabel('\theta [deg]')
title('angular distribution')
% simresult=simresult_array(11,9);
%%
indB_include=3;
Wronsk = cat(3,simresult.moments(length(simresult.delay0)+1:end,1:indB_include),simresult.moments_diff(:,1:indB_include,1:indB_include-1));
for indTT=1:size(Wronsk,1)
    Wronsk_t(indTT)=det(squeeze(Wronsk(indTT,:,:)));
end
figure;plot(simresult.delay(length(simresult.delay0)+1:end),Wronsk_t,'k')
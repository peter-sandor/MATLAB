in.FHCP=0; % laser peak intensity in [W/cm2]
in.I0=5e12; % laser peak intensity in [W/cm2]
in.TK=1; % rotational temperature in [K]
in.DELAY=0; % relative delay of orientation pulse with respect to alignment pulse
in.TAU=2; % Full cycle peak separation [ps]
in.TAUL=0.075; % laser pulse duration (intensity FWHM in [ps])
in.AD=31; % polarizability anisotropy [atomic units]
in.R0=0.33; % molecular dipole moment [atomic units]
in.B=0.2026; % rotational constant 'B' in [1/cm]
in.D=0*3.46e-8; % centrifugal distortion 'D' [1/cm]
in.EVEN=1; % abundance of even J states
in.ODD=1; % abundance of odd J states
in.JMAX=50; % max J level considered
in.THETANUM=200;
in.TMAX=90; % max timedelay in [ps]
in.DT=0.1;
in.ACCUR=1e-5;

simresult2=oalign3_matlab_v2(in);
%%
Jvec2=(0:in.JMAX+1);
amp0=squeeze(out.amps(end,:,:,:));
% amp0=out.amps_end;
amp1=amp0;
amp1(:,:,2:end)=sqrt(2)*squeeze(amp1(:,:,2:end));
amp1=squeeze(sum(amp1,2));
Prob1=squeeze(sum(abs(amp1).^2.*extend(out.AMP(1:size(amp0,1)),12),2));
sum(Prob1)
figure;hold on;
plot(Jvec2,out.AMP.*(2*Jvec2+1),'ko')
plot(Jvec2,Prob1,'ro')
plot(Jvec,Jpop,'bo')
% legend('oalign3 initial','oalign3 output','solve\_align code output')
%%
% load('I:\group\codes\Alignment_Codes\oalign3_linear\scan_I_and_T_75fs_OCS.mat')
figure;plot(out.TT*24e-6,out.TCOS2,'k')
% indT=9;indI=4;
% hold on;plot(simresult_array(indT,indI).delay,simresult_array(indT,indI).moments(:,3),'r')
% legend('MATLAB','FORTRAN')
% title(['T=' num2str(input_params(indT,indI).Trot) ' K; I=' num2str(input_params(indT,indI).laser_int*1e-12) ' TW/cm2'])
figure;imagescP(out.PROB)
%%
J=10;
M=5;
indJ=J+1;
M_offset=(size(simresult.amps,3)-1)/2+1;
figure;hold on;
plot(out.t_ip,abs(squeeze(out.amps(:,indJ,M+1,indJ))),'k')
plot(simresult.delay,abs(squeeze(simresult.amps(:,indJ,M_offset+M,indJ)))/abs(squeeze(simresult.amps(1,indJ,M_offset+M,indJ))),'r')
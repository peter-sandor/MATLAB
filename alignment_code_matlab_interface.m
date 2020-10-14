foldername='.\';
cd(foldername)

input_params.out_fname='out.dat';
input_params.outfull_fname='full_out.dat';
input_params.HCP_field=0.00001; % peak HCP field [kV/cm]
input_params.peaksep=2.0; % full cycle peak separation [ps]
input_params.laser_int=4.3e13; % laser peak intensity in [W/cm^2]
input_params.laser_fwhm=0.03; % laser pulse duration intensity FWHM in [ps]
input_params.laser_delay=0; % delay of laser after HCP in [ps]
input_params.polar_anisotropy=31; % polarizability anisotropy [atomic units]
input_params.dipole=0.33; % molecular dipole moment [atomic units]
input_params.rot_cnst=0.203; % rotational constant 'B' in [1/cm]
input_params.centrif=3.46e-8; % centrifugal distortion 'D' [1/cm]
input_params.abund_evenJ=1; % abundance of even J states
input_params.abund_oddJ=1; % abundance of odd J states
input_params.maxJ=90; % max J level considered
input_params.Trot=30; % rotational temperature in [K]
input_params.maxdelay=100; % max timedelay in [ps]
input_params.delaystep=0.1; % delay stepsize in [ps]
input_params.step_acc=1e-8; % single step accuracy

[simresult codeprint] = run_oalign3(input_params);

save simresult.mat simresult input_params;

fileID1 = fopen([foldername 'code_print.dat'],'w+');
fprintf(fileID1,'%s',codeprint);
fclose(fileID1);

%% Plot thermal J distribution
% rotational constants
B=input_params.rot_cnst; % [1/cm]
D=input_params.centrif; % [1/cm]
T=input_params.Trot;

J=0:input_params.maxJ;
h=6.626e-34; % [J*s]
c=3e8; % [m/s]
kB=1.38e-23; % [J/K]
E=h*c*100*(B*J.*(J+1)+D*J.^2.*(J+1).^2);
p=MB_distr(E,T);
coeffsJ=p/sum(p);

figure;
% subplot(121)
% plot(J,E/1.6e-19)
% subplot(122)
bar(J,coeffsJ)
xlabel('J')
ylabel('thermal population');
xlim([0 sqrt(kB*T/h/c/B)/3])
setfigP;
% 
%% plots the angular distribution output files from orient2d - the hcp orientation calculation
load simresult;
t0=min(simresult.delay);
tmax=max(simresult.delay);

figure;
subplot(231)
maxy=max(simresult.moments(:,1))*1.03;
miny=min(simresult.moments(:,1))/1.05;
plot(simresult.delay,simresult.moments(:,1),'LineWidth',2);
axis([0 tmax 0 maxy]);
   xlabel('Time (psec)','Fontsize',24);
   ylabel('<cos^2\theta> (2D)','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);
   axis ([t0 tmax miny maxy]); 
   axes('FontSize',20);
subplot(232)
maxy=max(simresult.moments(:,2))*1.03;
miny=min(simresult.moments(:,2))/1.05;
plot(simresult.delay,simresult.moments(:,2),'LineWidth',2);
axis([0 tmax 0 maxy]);
   xlabel('Time (psec)','Fontsize',24);
   ylabel('<cos^4\theta> (2D)','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);
   axis ([t0 tmax miny maxy]); 
   axes('FontSize',20);
subplot(233)
   surf(simresult.delay,simresult.theta,squeeze(simresult.data(:,:,1)).','FaceColor','interp','EdgeColor','none','AmbientStrength',.8)
   lighting phong;
%   camlight headlight;
   view(-90,85);
   axis manual;
   axis([t0 tmax 0 1 0 max(max(squeeze(simresult.data(:,:,1))))]);
   xlabel('Time (psec)','Fontsize',24);
   ylabel('\theta (\pi)','Fontsize',24);
   zlabel('Probability','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);

subplot(235)
maxy=max(simresult.moments(:,4))*1.03;
miny=min(simresult.moments(:,4))/1.05;
plot(simresult.delay,simresult.moments(:,4),'LineWidth',2)
   xlabel('Time (psec)','Fontsize',24);
   ylabel('<cos^4\theta> (3D)','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);
   axis ([t0 tmax miny maxy]);
subplot(234)
maxy=max(simresult.moments(:,3))*1.03;
miny=min(simresult.moments(:,3))/1.05;
plot(simresult.delay,simresult.moments(:,3),'LineWidth',2)
   xlabel('Time (psec)','Fontsize',24);
   ylabel('<cos^2\theta> (3D)','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);
   axis ([t0 tmax miny maxy]); 
subplot(236)
   surf(simresult.delay,simresult.theta,squeeze(simresult.data(:,:,2)).','FaceColor','interp','EdgeColor','none','AmbientStrength',.8)
   lighting phong;
%   camlight headlight;
   view(-90,85);
   axis manual;
   axis([t0 tmax 0 1 0 max(max(squeeze(simresult.data(:,:,1))))]);
   xlabel('Time (psec)','Fontsize',24);
   ylabel('\theta (\pi)','Fontsize',24);
   zlabel('Probability','Fontsize',24);
   set(gca,'LineWidth',2,'FontSize',20);
function [simresult codeprint] = run_oalign4(input_params)

Nlines=15;
foldername0=pwd;
foldername1='I:\group\codes\Alignment_Codes\oalign3_linear\';
cd(foldername1);

string_towrite{1}=input_params.out_fname;
string_towrite{2}=input_params.outfull_fname;
string_towrite{3}=num2str(input_params.HCP_field,'%.5f'); % peak HCP field [kV/cm]
string_towrite{4}=num2str(input_params.peaksep,'%.5f'); % full cycle peak separation [ps]
string_towrite{5}=num2str(input_params.laser_int,'%.5e'); % laser peak intensity in [W/cm^2]
string_towrite{6}=num2str(input_params.laser_fwhm); % laser pulse duration intensity FWHM in [ps]
string_towrite{7}=num2str(input_params.laser_delay,'%.5f'); % delay of laser after HCP in [ps]
string_towrite{8}=num2str(input_params.polar_anisotropy,'%.5f'); % polarizability anisotropy [atomic units]
string_towrite{9}=num2str(input_params.dipole); % molecular dipole moment [atomic units]
string_towrite{10}=[num2str(input_params.rot_cnst) ', ' num2str(input_params.centrif)]; % rotational constant 'B' and centrifugal distortion 'D' in [1/cm]
string_towrite{11}=[num2str(input_params.abund_evenJ,'%.1f') ', ' num2str(input_params.abund_oddJ,'%.1f')]; % abundance of even and odd J states
string_towrite{12}=num2str(input_params.maxJ); % max J level considered
string_towrite{13}=num2str(input_params.Trot); % rotational temperature in [K]
string_towrite{14}=[num2str(input_params.maxdelay,'%.5f') ',' num2str(input_params.delaystep)]; % max timedelay and delay stepsize in [ps]
string_towrite{15}=num2str(input_params.step_acc,'%.5e'); % single step accuracy

fileID1 = fopen([foldername1 'oalign-3.dat'],'w+');
for ind2=1:Nlines
    fprintf(fileID1,'%s\n',string_towrite{ind2});
end
fclose(fileID1);
%%
if exist([foldername1 input_params.out_fname])==2
    delete([foldername1 input_params.out_fname]);
end
if exist([foldername1 input_params.outfull_fname])==0
    fileID2=fopen([foldername1 input_params.outfull_fname],'w');
    fclose(fileID2);
end
tic;
[status,codeprint] = system([foldername1 'oalign-4.exe']); % run the calculation
simresult.runtime=toc;
%% Calculate moments <cos^2(theta)> and <cos^4(theta)> for 2D and 3D cases
simresult.pops=printconvert(codeprint);
data=load([foldername1 input_params.outfull_fname]);
norm=sum(data(:,2:end),2);
tstep=data(2,1)-data(1,1);
t0=data(1,1);
thstep = pi/(length(data(1,:))-1);
data(:,1)=[];
zmax=max(data);
zzmax=max(zmax);
datasize=size(data);
tnum=datasize(1,1);
thnum=datasize(1,2);
delay=(0:tnum-1);
delay=delay*tstep+t0;
tmax=max(delay);
theta=thstep*(0:thnum-1)/pi;
delay=delay';
xx=delay*(0*theta+1);
% theta=(0*delay+1)*theta;
theta_ext=permute(extend(theta,length(delay)),[2 1]);

cos_sq_2D=sum(data.*cos(pi*theta_ext).^2,2)./norm;   % cos^2 2D
cos_quad_2D=sum(data.*cos(pi*theta_ext).^4,2)./norm; %cos^4 2D
cos_sq_3D=sum(data.*sin(pi*theta_ext).*cos(pi*theta_ext).^2,2)*thstep;
cos_quad_3D=sum(data.*sin(pi*theta_ext).*cos(pi*theta_ext).^4,2)*thstep;

simresult.input=input_params;
simresult.theta=map2colvec(theta);
simresult.delay=map2colvec(delay);
% simresult.data=cat(3,data,data2);
simresult.data=data;
% simresult.moments=[map2colvec(cos_sq_2D) map2colvec(cos_quad_2D) map2colvec(cos_sq_3D) map2colvec(cos_quad_3D)];
simresult.moments=[map2colvec(cos_sq_3D) map2colvec(cos_quad_3D)];
cd(foldername0);
end
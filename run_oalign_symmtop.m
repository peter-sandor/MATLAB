function [simresult codeprint] = run_oalign_symmtop(input_params)

Nlines=16;
foldername0=pwd;
foldername1='I:\group\codes\Alignment_Codes\align_symmetric_top\';
cd(foldername1);
fileID1 = fopen([foldername1 'Align.dat'],'w+');

string_towrite{1}=input_params.out_fname;
string_towrite{2}=input_params.outfull_fname;
string_towrite{3}=num2str(input_params.HCP_field,'%.5f'); % peak HCP field [kV/cm]
string_towrite{4}=num2str(input_params.peaksep,'%.5f'); % full cycle peak separation [ps]
string_towrite{5}=num2str(input_params.laser_int2,'%.5e'); % laser peak intensity in [W/cm^2]
string_towrite{6}=num2str(input_params.laser_fwhm2); % laser pulse duration intensity FWHM in [ps]
string_towrite{7}=num2str(input_params.t0,'%.3f'); % delay of laser after HCP in [ps]
string_towrite{8}=num2str(input_params.polar_antr,'%.3f'); % polarizability anisotropy [atomic units]
string_towrite{9}=num2str(input_params.dipole); % molecular dipole moment [atomic units]
string_towrite{10}=[num2str(input_params.rot_cnst_B,'%.5f') ', ' num2str(input_params.centrif,'%.3e')]; % rotational constant 'B' and centrifugal distortion 'D' in [1/cm]
string_towrite{11}=[num2str(input_params.abund_evenJ,'%.1f') ', ' num2str(input_params.abund_oddJ,'%.1f')]; % abundance of even and odd J states
string_towrite{12}=num2str(input_params.maxJ); % max J level considered
string_towrite{13}=num2str(input_params.Trot); % rotational temperature in [K]
string_towrite{14}=[num2str(input_params.maxdelay,'%.3f') ',' num2str(input_params.delaystep,'%.3f')]; % max timedelay and delay stepsize in [ps]
string_towrite{15}=num2str(input_params.step_acc,'%.3e'); % single step accuracy
string_towrite{16}=num2str(input_params.rot_cnst_A,'%.5f'); % rotational constant 'A' in [1/cm]

% fileID2=fopen('oalign-3.dat','r+');
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
if exist([foldername1 input_params.jamp_fname])==0
    fileID2=fopen([foldername1 input_params.jamp_fname],'w');
    fclose(fileID2);
end
code_to_run='oalign_symmtop5.exe';
% code_to_run='oalign_symmtop_test.exe';
[status,codeprint] = system([foldername1 code_to_run]); % run the calculation
% status=[];
% codeprint=[];
%%
probs0=dlmread('jamp.dat');
simresult.probs_J=probs0(:,2);
if strcmp(code_to_run,'oalign_symmtop_test.exe')
    amps_re=dlmread('amplitude-re.dat');
    amps_im=dlmread('amplitude-im.dat');
    Mvals=unique(amps_re(:,1));
    Kvals=Mvals;
    maxM=max(Mvals);
    maxK=maxM;
    amps=zeros([size(amps_re,2)-3 2*maxM+1 2*maxK+1 maxM+1]);
    ind_offset=maxM+1;
    for ind1=1:size(amps_re,1)
        amps(:,amps_re(ind1,1)+ind_offset,amps_re(ind1,2)+ind_offset,amps_re(ind1,3)+1) = map2colvec(amps_re(ind1,4:end)) + 1i*map2colvec(amps_im(ind1,4:end));
    end
    simresult.amps=amps;
end
%% Calculate moments <cos^2(theta)> and <cos^4(theta)> for 2D and 3D cases
simresult.input=input_params;
simresult.codeprint=codeprint;
if exist([foldername1 input_params.out_fname])==2
    fmom=load([foldername1 input_params.out_fname]);
else
    fmom=zeros([1 3]);
end
% delay=fmom(:,1);
data=load([foldername1 input_params.outfull_fname]);
delay=data(:,1);
data(:,1)=[];
norm=sum(data,2);

thstep = 1/(size(data,2)-1);
[Ndelay,Ntheta]=size(data);
theta=thstep*(0:Ntheta-1);

simresult.input=input_params;
simresult.theta=map2colvec(theta);
simresult.delay=map2colvec(delay);
simresult.tcos1=fmom(:,2);
simresult.tcos2=fmom(:,3);
simresult.data=data;
cd(foldername0);
end
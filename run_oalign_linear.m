function [simresult codeprint] = run_oalign_linear(input_params)

foldername0=pwd;
foldername1='I:\group\codes\Alignment_Codes\oalign3_linear\';
cd(foldername1);

string_towrite{1}=input_params.out_fname;
string_towrite{2}=input_params.outfull_fname;
string_towrite{3}=num2str(input_params.HCP_field,'%.5f'); % peak HCP field [kV/cm]
string_towrite{4}=num2str(input_params.peaksep,'%.5f'); % full cycle peak separation [ps]
string_towrite{5}=num2str(input_params.laser_int2,'%.5e'); % laser peak intensity in [W/cm^2]
string_towrite{6}=num2str(input_params.laser_fwhm2); % laser pulse duration intensity FWHM in [ps]
string_towrite{7}=num2str(input_params.t0,'%.5f'); % delay of laser after HCP in [ps]
string_towrite{8}=num2str(input_params.polar_antr,'%.5f'); % polarizability anisotropy [atomic units]
string_towrite{9}=num2str(input_params.dipole); % molecular dipole moment [atomic units]
string_towrite{10}=[num2str(input_params.rot_const_B) ', ' num2str(input_params.centrif)]; % rotational constant 'B' and centrifugal distortion 'D' in [1/cm]
string_towrite{11}=[num2str(input_params.abund_evenJ,'%.1f') ', ' num2str(input_params.abund_oddJ,'%.1f')]; % abundance of even and odd J states
string_towrite{12}=num2str(input_params.maxJ); % max J level considered
string_towrite{13}=num2str(input_params.Trot); % rotational temperature in [K]
string_towrite{14}=[num2str(input_params.maxdelay,'%.5f') ',' num2str(input_params.delaystep)]; % max timedelay and delay stepsize in [ps]
string_towrite{15}=num2str(input_params.step_acc,'%.5e'); % single step accuracy
% string_towrite{16}=num2str(input_params.rot_cnst_A,'%.5f'); % rotational constant 'A' in [1/cm]
Nlines=length(string_towrite);

% fileID1 = fopen([foldername1 'oalign-3.dat'],'w+');
fileID1 = fopen([foldername1 'oalign-5.dat'],'w+');
for ind2=1:Nlines
    fprintf(fileID1,'%s\n',string_towrite{ind2});
end
fclose(fileID1);

% fileID3 = fopen([foldername1 'jamp.dat'],'w+');
% fprintf(fileID1,'%s','');
% fclose(fileID3);
%%
if exist([foldername1 input_params.out_fname])==2
    delete([foldername1 input_params.out_fname]);
end
if exist([foldername1 input_params.outfull_fname])==0
    fileID2=fopen([foldername1 input_params.outfull_fname],'w');
    fclose(fileID2);
end
% code_to_run='oalign_linear7_test.exe';
code_to_run='oalign_linear7.exe';
tic;
[status,codeprint] = system([foldername1 code_to_run]); % run the calculation
simresult.runtime=toc;
simresult.pops=printconvert(codeprint);
simresult.codeprint=codeprint;

if strcmp(code_to_run,'oalign_linear7_test.exe')
    amps_re=dlmread('amplitude-re.dat');
    amps_im=dlmread('amplitude-im.dat');
    Mvals=unique(amps_re(:,1));
    maxM=max(Mvals);
    amps=zeros([size(amps_re,2)-2 maxM+1 maxM+1]);
    for ind1=1:size(amps_re,1)
        amps(:,amps_re(ind1,1)+1,amps_re(ind1,2)+1) = map2colvec(amps_re(ind1,3:end)) + 1i*map2colvec(amps_im(ind1,3:end));
    end
    simresult.amps=amps;
end
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
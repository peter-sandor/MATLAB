function FROG2ifa(varargin)

if nargin==1
    [folder_name,FileName]=strip_path(varargin{1});
else
    folder_name=[pwd '\'];
    FileName='FROG.txt';
%     [FileName,folder_name,FilterIndex] = uigetfile('*.txt');
end
%% load files
FROG=load([folder_name FileName]);
delays=load([folder_name 'delays.txt']);
wavelengths=load([folder_name 'wavelengths.txt']);
[delays_con,FROG_con,ind1]=consolidator(map2colvec(delays),FROG,'mean',0.0001);
delays=delays_con;
FROG=FROG_con;
%% determine timezero, center wavelength and stepsize
Ndelays=length(delays);
Nwl=length(wavelengths);
N=256;

T0=sum(squeeze(sum(FROG,2)).*map2colvec(delays))/sum(sum(FROG));
delays_fs=(delays-T0)*1e6/299.792*2;
%% subtract background and threshold
Nbkg1=round(sqrt((Ndelays*0.05)^2+(Nwl*0.05)^2));
Nbkg2=Nbkg1;
minval=min(min(FROG));
maxval=max(max(FROG));
bkg_vector=mean(FROG(1:Nbkg2,:),1);
bkg_mean=mean(mean(FROG(1:Nbkg1,1:Nbkg1)));
bkg_max=max(max(FROG(1:Nbkg1,1:Nbkg1)));
bkg_min=min(min(FROG(1:Nbkg1,1:Nbkg1)));
FROG=FROG-permute(extend(bkg_vector,Ndelays),[2 1]);
thrs=1.0*(bkg_max-bkg_mean);
FROG(FROG<=thrs)=0;

delays_ip=min(delays_fs):(max(delays_fs)-min(delays_fs))/(N-1):max(delays_fs);
wl_ip=min(wavelengths):(max(wavelengths)-min(wavelengths))/(N-1):max(wavelengths);
delay_step=abs(mean(diff(delays_ip)));
wl_step=mean(diff(wl_ip));
wl_center=wl_ip(ceil(N/2));
[TAU_ip,WL_ip]=meshgrid(delays_ip,wl_ip);
[TAU,WL]=meshgrid(delays_fs,wavelengths);
FROG_ip = interp2(TAU,WL,FROG.',TAU_ip,WL_ip);
%% write to .ifa
ifa_fullpath=[folder_name 'FROG.ifa'];
fileID=fopen(ifa_fullpath,'w');
fprintf(fileID,'%d\n',N);
fprintf(fileID,'%d\n',N);
fprintf(fileID,'%2.7E\n',delay_step);
fprintf(fileID,'%3.7E\n',wl_step);
fprintf(fileID,'%3.7E\n',wl_center);
fclose(fileID);
dlmwrite(ifa_fullpath,FROG_ip,'precision','%1.3E','delimiter','\t','newline','pc','-append');
end
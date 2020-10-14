function filter_FROG(varargin)
if nargin==1
    folder_name=varargin{1};
else
    folder_name = [uigetdir '\'];
end
%% load files
FROG=load([folder_name 'FROG.txt']);
delays=load([folder_name 'delays.txt']);
wavelengths=load([folder_name 'wavelengths.txt']);
%% sutract background and threshold
Nbkg1=50;
Nbkg2=10;
Ndelays=length(delays);
Nwl=length(wavelengths);

bkg_vector=mean(FROG(1:Nbkg2,:),1);
bkg_mean=mean(mean(FROG(1:Nbkg1,1:Nbkg1)));
bkg_max=max(max(FROG(1:Nbkg1,1:Nbkg1)));
bkg_min=min(min(FROG(1:Nbkg1,1:Nbkg1)));
FROG=FROG-permute(extend(bkg_vector,Ndelays),[2 1]);
thrs=ceil(bkg_max-bkg_mean);
FROG(FROG<=thrs)=0;

T0=sum(squeeze(sum(FROG,2)).*map2colvec(delays))/sum(sum(FROG));
delays_fs=(delays-T0)*1e6/299.792*2;
%% compute FFT and filter
freq_dels=FourierAxis(delays_fs);
FROG_FFT=fft(FROG,[],1);
% figure;imagescP(wavelengths,fftshift(freq_dels),fftshift(abs(FROG_FFT),1))
fhndl1=figure;imagescP(abs(FROG_FFT))
temp=round(ginput(1));
index1=temp(2);
maskwidth=0.5;
maskfun=1./(exp((map2colvec(1:Ndelays) - index1)/maskwidth)+1);
mask=extend(maskfun+flipdim(maskfun,1),Nwl);
close(fhndl1)
FROG_FFT_filtered=FROG_FFT.*mask;
FROG_filtered=ifft(FROG_FFT_filtered,[],1);
figure;imagescP(delays_fs,wavelengths,abs(FROG_filtered).')
dlmwrite([folder_name 'FROG_filtered.txt'],FROG_filtered,'\t');
end
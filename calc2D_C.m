function varargout = calc2D_C(varargin)

% Processing
% This code handles data taken at single pump-probe delays, read in with 'dataload_pscan.m'
% the dimensions of 'data' should be the following: [nrphases nrppdelays nrchannels Nwl]
% Subtract background, interpolate to have equidistant spacing in frequency space, divide probe channel by ref boxcar, take difference of
% phase cycles, then multiply signal back with probe spectra averaged over all pump-pump delays
% procpath='.\';
c=300; % [nm/fs]
if nargin==0
    load([procpath 'data']);
    procpath='.\';
elseif nargin==1
    data=varargin{1};
    procpath='.\';
elseif nargin==2
    data=varargin{1};
    procpath=varargin{2};
else 'Incorrect number of arguments.'
end

load([procpath 'Np']);
load([procpath 'axes']);
if exist('frame.mat')==2
    load([procpath 'frame']);
    if length(frame)==0
        answer=inputdlg('Enter rotating frame frequency in PHz:');
        frame=str2num(answer{1});
    end
else answer=inputdlg('Enter rotating frame frequency in PHz:');
    frame=str2num(answer{1});
end

if length(size(data))==5
    data=squeeze(mean(data,3));
end
data=flipdim(data,2);
ppdelay_max=max(axes{2}); % [fs]
nrchannels=2;
nrphases=size(data,1);
nrppdelays=size(data,2);
ZFF=4;
N1=ZFF*nrppdelays;
N2=128;
frq1_step=1/ppdelay_max/ZFF;
if mod(N1,2)==0
    frq1_min=-N1/2*frq1_step;
    frq1_max=(N1/2-1)*frq1_step;
elseif mod(N1,2)==1
    frq1_min=-1/(2*delay_step);
    frq1_max=1/(2*delay_step);
end
frq1=map2colvec(0:frq1_step:frq1_max)+frame;

Nwl=size(data,4);
if exist('..\w1')==2
    w1=load('..\w1');
    w2=load('..\w2');
else w1=load('..\w1.dat');
    w2=load('..\w2.dat');
end
wl_max=min([max(w1),max(w2)]);
wl_min=max([min(w1),min(w2)]);

freqch1=map2colvec(c./w1);
freqch2=map2colvec(c./w2);
frq3=map2colvec(c/wl_max:frq1_step:c/wl_min);
N3=size(frq3,1);
%%
bkg_begin=1;
bkg_end=50;
box1=1;
box2=N3;

data_ip=zeros([nrphases nrppdelays nrchannels N3]);
for ind3=1:nrppdelays
	for ind2=1:nrphases
        data(ind2,ind3,1,:)=data(ind2,ind3,1,:)-mean(data(ind2,ind3,1,bkg_begin:bkg_end),4);
        data_ip(ind2,ind3,1,:)=squeeze(interp1(freqch1,squeeze(data(ind2,ind3,1,:)),frq3,'linear'));
        data(ind2,ind3,2,:)=data(ind2,ind3,2,:)-mean(data(ind2,ind3,2,bkg_begin:bkg_end),4);
        data_ip(ind2,ind3,2,:)=squeeze(interp1(freqch2,squeeze(data(ind2,ind3,2,:)),frq3,'linear'));
        data_div(ind2,ind3,:)=squeeze(data_ip(ind2,ind3,1,:)./(abs(data_ip(ind2,ind3,2,:))+1));
    end
    data_diff(ind3,:)=data_div(2,ind3,:)-data_div(1,ind3,:);
end
% data_diff=data_diff2.*(ones([nrppdelays 1])*squeeze(mean(data_ip(1,1,:,:),3)).');
% remove DC component
data_diff=data_diff-mean(mean(data_diff));
data_diff=flipdim(data_diff,1);
clear data_diff2
%% Fourier transform along pump-pump axis
% zeropadding for higher spectral resolution (also helps avoiding "picket fence"-effect)
marg1=squeeze(mean(data_diff,2));
data_diff2=[data_diff; zeros([N1-nrppdelays N3])];
data_2D_full=fftshift(fft(data_diff2,[],1),1);
% The way to set up the frequency axis and to select the data corresponding to
% positive frequencies depends on whether the array that was transformed contains even or odd
% number of elements
if mod(N1,2)==0
    data2D=data_2D_full(N1/2+1:N1,:);
elseif mod(N1,2)==1
    data2D=data_2D_full((N1+1)/2+1:N1,:);    
end
varargout{1}=data2D;
varargout{2}=frq3;
varargout{3}=frq1;
varargout{4}=marg1;
end
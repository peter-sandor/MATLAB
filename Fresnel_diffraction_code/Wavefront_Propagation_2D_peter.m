%% Square Aperture

lambda=1; % [micron]
k=2*pi/lambda;
padfactor=1;
interpfactor=1;
phasefactor=0;
f=1000*1000;
n=1.5;

x0range=5000;
Nx0=1024;
x0step=x0range/(Nx0-1);
x0=-x0range/2:x0step:x0range/2;
[X0,Y0]=meshgrid(x0,x0);

sidelength=1000;
U0=zeros(size(X0));
U0(logical((abs(X0)<sidelength/2).*(abs(Y0)<sidelength/2)))=1;

%x1range^2/5/z^3/lambda; % this quantity has to be << 1 for Fresnel-approximation to be valid everywhere on the screen
z_min=5*(x0range^2/5/lambda)^(1/3); % alternatively, the minimum 'z' can be defined based on the above
% << is approximated by 1% (x<<y <--> x<0.01*y).
z=[1:2:8]*16*z_min; % Nx0=1024, lambda=1 --> z>=16*min_z

%% Circular Aperture
lambda=1; % [micron]
k=2*pi/lambda;
padfactor=1;
interpfactor=1;
sigmaF=20;
phasefactor=0;
f=1000*1000;
n=1.5;

x0range=5000;
Nx0=1024;
x0step=x0range/(Nx0-1);
x0=-x0range/2:x0step:x0range/2;
[X0,Y0]=meshgrid(x0,x0);

%x1range^2/5/z^3/lambda; % this quantity has to be << 1 for Fresnel-approximation to be valid everywhere on the screen
z_min=5*(x0range^2/5/lambda)^(1/3); % alternatively, the minimum 'z' can be defined based on the above
% << is approximated by 1% (x<<y <--> x<0.01*y).
z=[1:2:8]*16*z_min; % Nx1=1024, lambda=1 --> z>=16*min_z

circ_radius_index=1*round(size(X0,1)*1/8);
% U0=zeros(size(X1));
% U0=exp(-(X1.^2+Y1.^2)/2/sigmaF^2).*double(imagecirc(size(X1),round(size(X1)/2),[0*round(size(X1,1)*1/8) circ_radius_index]));
U0=1.*double(imagecirc(size(X0),round(size(X0)/2),[0*round(size(X0,1)*1/8) circ_radius_index]));
% U_analytic=besselj(1,pi*2*circ_radius_index*x0step*sqrt(X0.^2+Y0.^2)/lambda/z)./(pi*2*circ_radius_index*x0step*sqrt(X0.^2+Y0.^2)/lambda/max(z)); % analytic solution for circular aperture in the Fraunhofer limit

%% Gaussian beam
lambda=1; % [micron]
k=2*pi/lambda;
sigma=2000;
padfactor=2;
interpfactor=16;
phasefactor=1;
f=1000*1000;
n=1.5;

x0range=20*1000;
Nx0=128;
x0step=x0range/(Nx0-1);
x0=-x0range/2:x0step:x0range/2;
[X0,Y0]=meshgrid(x0,x0);
U0=exp(-(X0.^2+Y0.^2)/sigma^2);

%% Amplifier beam profile

% uiopen('M:\2015_04_10_E\amp_beam_image.fig',1)
% [x,y,intensity_profile]=getdata;
load amplifier_profile.mat;
lambda=1; % [micron]
k=2*pi/lambda;  
padfactor=1;
interpfactor=16;
phasefactor=1;
f=1000*1000;
n=1.5;

x0=x*1000; % [x]=mm, [x0]=micron
Nx0=length(x0);
x0range=max(x0)-min(x0);
x0step=x0range/(Nx0-1);
[X0,Y0]=meshgrid(x0,x0);
mask=double(abs(Y0)>=x0range*(-1));
U0=sqrt(intensity_profile).*mask; % we need the field
%% Preprocess: interpolate and zeropad
% First interpolate input
x1step=x0range/(Nx0*interpfactor-1);
x1=-x0range/2:x1step:x0range/2;
y1=x1;
[X1,Y1]=meshgrid(x1,x1);
U1=interp2(X0,Y0,U0,X1,Y1);
Nx1=size(U1,1);

% Next, zeropad
x2step=padfactor*x0range/(padfactor*Nx1-1);
x2=-padfactor*x0range/2:x2step:padfactor*x0range/2;
Nx2=length(x2);
[X2,Y2]=meshgrid(x2,x2);
U2=zeros(padfactor*[Nx1 Nx1]);
ind_insert=round((padfactor-1)/2*Nx1(1)+1:1:(padfactor+1)/2*Nx1);
U2(ind_insert,ind_insert)=U1;

% write phase for a plano-convex lens:

R=f*(n-1);
phase=phasefactor*k*(n-1)*(sqrt(R^2-(X2.^2+Y2.^2))-1*72e3*ones(size(X2)))+pi*double(Y2>0); % take positive sign
U_object=U2.*exp(1i*phase);
z=f;%*(0.8:0.1:1.2);
%% Do the calculation
Nz=length(z);
U=zeros([Nx2 Nx2 Nz]);
for ind1=1:Nz
    U_image = Fresnel_propagate_CT(U_object,lambda,z(ind1),x2);
%     [U2,b,c]=Fresnel_propagate_basic(U0,lambda,z(ind1),x1,x2);
    U(:,:,ind1) = U_image;
    disp([num2str(ind1) '/' num2str(Nz)])
end

figure(1);
subplot(1,Nz+1,1)
imagesc(x2/1000,x2/1000,abs(U_object).^2)
% xlim([-10 10]);
% ylim([-10 10]);
xlabel('x [mm]')
ylabel('y [mm]')
title('|U_{object}|^2')
for ind1=1:Nz
    subplot(1,Nz+1,ind1+1)
    imagesc(x2/1000,x2/1000,abs(U(:,:,ind1)).^2)
%     xlim([-10 10]);
%     ylim([-10 10]);
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['|U_{image}|^2; z=' num2str(z(ind1)/1000) ' mm'])
end
% subplot(133)
% imagesc(x1/1000,y1/1000,abs(U_analytic).^2)
% xlabel('x [mm]')
% ylabel('y [mm]')
% title(['analytic solution; z=' num2str(z(index)/1000) ' mm'])
%%
figure;plot(x2,(squeeze(sum(abs(U(:,round(Nx2/2),1)),2))/max(squeeze(sum(abs(U(:,round(Nx2/2),1)),2)))).^2,'r');
hold on;plot(x2,(squeeze(sum(abs(U_object(:,round(Nx2/2),1)),2))/max(squeeze(sum(abs(U_object(:,round(Nx2/2),1)),2)))).^2,'k');
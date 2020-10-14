% test curve recognition
pic=monopic(imread('1D_potentials_from_abstract_cropped.bmp'));
pic=255-pic;
[x,y]=recognize_curve(pic);
xax=-0.8:1.6/293:0.8;
yax=7.25:1.75/146:8.5;
figure;plot(xax,yax(round(y)))
%%
pic=color2gray(imread('U:\People\S.Peter\projects\IFROG\BG23_cropped.bmp'));
pic=max(max(pic))-pic;
[x,y]=recognize_curve(double(pic),215);
x_min0 = 200; % [nm]
x_max0 = 1100; % [nm]
Nx0 = size(pic,2);
Ny0 = size(pic,1);
Nx = length(x);
xax=x_min0:(x_max0 - x_min0)/(Nx0-1):x_max0;

y_min = 0;
y_max = 100;
yax=y_min:(y_max-y_min)/(Ny0-1):y_max;
figure;plot(xax(round(x)),yax(round(y)))

x_ip = 200:900/(511):1100;
T_BG23 = interp1(xax(round(x)),yax(round(y)),x_ip,'linear');
T_BG23(isnan(T_BG23)) = 1e-2;
figure; hold on;
plot(xax(round(x)),yax(round(y)),'k')
plot(x_ip,T_BG23,'r')
%%
pic=color2gray(imread('U:\People\S.Peter\projects\IFROG\Layertec_103366_Z0509040_center_net_cropped_redrawn.bmp'));
pic=max(max(pic))-pic;
[x,y]=recognize_curve(double(pic),50);
x_min0 = 620; % [nm]
x_max0 = 980; % [nm]
Nx0 = size(pic,2);
Ny0 = size(pic,1);
Nx = length(x);
xax=x_min0:(x_max0 - x_min0)/(Nx0-1):x_max0;

y_min = -200;
y_max = 200;
yax=y_min:(y_max-y_min)/(Ny0-1):y_max;
figure;plot(xax(round(x)),yax(round(y)))

Nip = 256;
x_min1 = min(xax(round(x)));
x_max1 = max(xax(round(x)));
x_ip = x_min1:(x_max1-x_min1)/(Nip-1):x_max1;
spectrum = interp1(xax(round(x)),yax(round(y)),x_ip,'linear');
x_ip(isnan(spectrum)) = [];
spectrum(isnan(spectrum)) = [];
figure; hold on;
plot(xax(round(x)),yax(round(y)),'k')
plot(x_ip,spectrum,'r')
dlmwrite('U:\People\S.Peter\projects\IFROG\Layertec_103366_Z0509040_GDD_vs_lambda.txt',[map2colvec(x_ip) map2colvec(spectrum)],'\t');
%%
pic=color2gray(imread('U:\People\S.Peter\projects\IFROG\USB4000_CCD_spectral_response_cleaned.png'));
pic=max(max(pic))-pic;
[x,y]=recognize_curve(double(pic),50);
x_min0 = 400; % [nm]
x_max0 = 1200; % [nm]
Nx0 = size(pic,2);
Ny0 = size(pic,1);
Nx = length(x);
xax=x_min0:(x_max0 - x_min0)/(Nx0-1):x_max0;

y_min = 0;
y_max = 1;
yax=y_min:(y_max-y_min)/(Ny0-1):y_max;
figure;plot(xax(round(x)),yax(round(y)))

Nip = 256;
x_min1 = min(xax(round(x)));
x_max1 = max(xax(round(x)));
x_ip = x_min1:(x_max1-x_min1)/(Nip-1):x_max1;
spectrum = interp1(xax(round(x)),yax(round(y)),x_ip,'linear');
x_ip(isnan(spectrum)) = [];
spectrum(isnan(spectrum)) = [];
figure; hold on;
plot(xax(round(x)),yax(round(y)),'k')
plot(x_ip,spectrum,'r')
dlmwrite('U:\People\S.Peter\projects\IFROG\USB4000_CCD_spectral_response.txt',[map2colvec(x_ip) map2colvec(spectrum)],'\t');
%%
pic=color2gray(imread('U:\People\S.Peter\projects\IFROG\Venteon_spectrum_spec.bmp'));
pic=max(max(pic))-pic;
[x,y]=recognize_curve(double(pic),50);
x_min0 = 550; % [nm]
x_max0 = 1250; % [nm]
Nx0 = size(pic,2);
Ny0 = size(pic,1);
Nx = length(x);
xax=x_min0:(x_max0 - x_min0)/(Nx0-1):x_max0;

y_min = 0;
y_max = 1;
yax=y_min:(y_max-y_min)/(Ny0-1):y_max;
figure;plot(xax(round(x)),yax(round(y)))

Nip = 256;
x_min1 = min(xax(round(x)));
x_max1 = max(xax(round(x)));
x_ip = x_min1:(x_max1-x_min1)/(Nip-1):x_max1;
spectrum = interp1(xax(round(x)),yax(round(y)),x_ip,'linear');
x_ip(isnan(spectrum)) = [];
spectrum(isnan(spectrum)) = [];
figure; hold on;
plot(xax(round(x)),yax(round(y)),'k')
plot(x_ip,spectrum,'r')
dlmwrite('U:\People\S.Peter\projects\IFROG\Venteon_spectrum_from_specs.txt',[map2colvec(x_ip) map2colvec(spectrum)],'\t');
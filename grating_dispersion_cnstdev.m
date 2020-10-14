% calculation of linear dispersion caused by a f = 250 mm monochromator equipped with a d = 1200
% groove/mm grating
% reference:
% http://www.tau.ac.il/~phchlab/experiments/hydrogen/diffraction_gratings.h
% tm
% Newport: Diffraction Grating Handbook

m = +2;
G = 1800; % [groove/mm]
d = 1/G; % [mm/groove]
d_SI = 1/(G*1e3); % [m/groove]
N=100;
% lambda_min=2.3e-7;
% lambda_max=2.9e-7;
lambda_min=2.3e-7;
lambda_max=2.9e-7;
delta_lambda = (lambda_max-lambda_min)/(N-1);
lambda = lambda_min:delta_lambda:lambda_max; % diffracted wavelengths in [m]

f = 75; % [mm]
Kappa = 15;
Kappa_rad = Kappa/180*pi;

% SONY ILX511: linear array of 2048 pixels, 14 um/pixel, 28.6720 mm total length of array

Phi = asin(m*lambda/(2*d_SI*cos(Kappa_rad))) *180/pi;
Alpha = Kappa+Phi;
Beta = Phi-Kappa;
beta_g = asin(m*lambda./d_SI-sin(Alpha/180*pi));
% % diffraction angle in Littrow configuration
% beta_L = asin(m*lambda/(2*d_SI));

% angular dispersion [d(beta)/d(lambda)]
ang_disp = (sin(Alpha/180*pi)+sin(Beta/180*pi))./(lambda.*cos(Beta/180*pi)); % [rad/m]

% ang_disp_L = m./(2*d_SI*cos(beta_L));
% ang_disp2_L = m/(2*d_SI)*1./sqrt(1-(m*lambda./(2*d_SI)).^2);

% reciprocal linear dispersion
rld = d*cos(Beta*pi/180)/(abs(m)*f);
rld_g = d*cos(beta_g)/(abs(m)*f);

% figure;
% [AX, H1, H2] = plotyy (lambda*1e9, rld*1e6, lambda*1e9, ang_disp*180*1e-9/pi);
% set (get(AX(1), 'Ylabel'), 'String', 'Reciprocal linear dispersion [nm/mm]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'magenta');
% set (get(AX(2), 'Ylabel'), 'String', 'Angular dispersion [degrees/nm]', 'fontname', 'verdana', 'fontsize', 12, 'color', 'blue');
% set (AX(1), 'Ycolor', 'magenta');
% set (AX(2), 'Ycolor', 'blue');
% set(H1, 'Linestyle', '-', 'Linewidth', 1, 'Color', 'magenta');
% set(H2, 'Linestyle', '-', 'Linewidth', 1, 'Color', 'blue');
% % set(AX,'xlim',[nu_lim1, nu_lim2]);
% %xlim ([nu_lim1 nu_lim2]);
% xlabel('Wavelength [nm]', 'fontname', 'verdana', 'fontsize', 12);
% title([num2str(G),' groove/mm;',' f = ' num2str(f), ' mm'], 'Fontsize', 14)

figure;
[AX, H1, H2] = plotyy (lambda*1e9, rld*1e6, lambda*1e9, ang_disp*180*1e-9/pi);
set (get(AX(1), 'Ylabel'), 'String', 'Reciprocal linear dispersion [nm/mm]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'magenta');
set (get(AX(2), 'Ylabel'), 'String', 'Angular dispersion [degrees/nm]', 'fontname', 'verdana', 'fontsize', 12, 'color', 'blue');
set (AX(1), 'Ycolor', 'magenta','Ylim',[min(rld)*1e6,max(rld)*1e6],'Ytick',[min(rld)*1e6:(max(rld)-min(rld))*1e6/4:max(rld)*1e6]);
set (AX(2), 'Ycolor', 'blue','Ylim',[min(ang_disp)*180*1e-9/pi,max(ang_disp)*180*1e-9/pi],'Ytick',[min(ang_disp)*180*1e-9/pi:(max(ang_disp)-min(ang_disp))*180*1e-9/pi/4:max(ang_disp)*180*1e-9/pi]);
set(H1, 'Linestyle', '-', 'Linewidth', 1, 'Color', 'magenta');
set(H2, 'Linestyle', '-', 'Linewidth', 1, 'Color', 'blue');
% set(AX,'xlim',[nu_lim1, nu_lim2]);
%xlim ([nu_lim1 nu_lim2]);
xlabel(AX(1),'diffracted wavelength [nm]', 'fontname', 'verdana', 'fontsize', 12);
text(lambda(1)*1e9,(.1*(max(rld)-min(rld))+min(rld))*1e6,[num2str(G),' groove/mm;',' f = ' num2str(f), ' mm;',' spectral order: ',num2str(m)], 'Fontsize', 14)

AXpos=get(AX,'Position');
ax2 = axes('Position',AXpos{1},'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor',[.8,.8,.8],'Visible','on','Ytick',[]);
hl2 = line(lambda*1e9*abs(m), rld*1e6,'Color','g','Parent',ax2,'LineStyle','none');
xlabel(ax2,'measured wavelength [nm]', 'fontname', 'verdana', 'fontsize', 12);
% If the ruling of the grating is 1200/mm, and the focal length is f = 240
% mm, an individual stage of the monochromator produces 3.4 nm/mm
% reciprocal linear dispersion. Two identical cascaded stages in the additive mode
% (CVI DK242) then produces 1(1/3.4 + 1/3.4) = 3.4/2 = 1.7 nm/mm.
% Provided the entrance slit width is equal to the detector pixel width and is equal to 11 micron (assuming that the fill factor is 100%),
% and the imaging is 1:1, one can achieve a spectral resolution of 0.2
% Angstrom. If the slit width is greater than the pixel width by a factor
% of two, then the spectral resolution is halved (0.4 Angstrom)
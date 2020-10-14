function varargout = prism_compressor(lambda0,apex_angle,theta0,a,t,s)

% calculate second and third order phase for a pair of FS prisms at any input angle, any prism apex angle, for a single pass.
% From Diels, Rudolph: Ultrashort Laser Pulse Phenomena

% apex_angle = prism apex angle in [deg]
% a,t,s = distances defined as in the book in [mm]
% omega in [rad/fs]
c=300e-6; % [mm/fs]

N=5;
omega0=c./lambda0*1e3*2*pi; % [rad/fs]
delta_omega=omega0*1e-3;
omega=map2colvec((omega0-(N-1)/2*delta_omega):delta_omega:(omega0+(N-1)/2*delta_omega));
omega_p=map2colvec((omega0-(N-2)/2*delta_omega):delta_omega:(omega0+(N-2)/2*delta_omega));
omega_pp=map2colvec((omega0-(N-3)/2*delta_omega):delta_omega:(omega0+(N-3)/2*delta_omega));
omega_ppp=map2colvec((omega0-(N-4)/2*delta_omega):delta_omega:(omega0+(N-4)/2*delta_omega));

refr_omega=sellmeier_FS(c./(omega/2/pi).*1000); % input to sellmeier_FS: [um]
refr_omega_p=diff(refr_omega)/delta_omega; % first derivative of the refractive index with respect to omega
refr_omega_pp=diff(refr_omega_p)/delta_omega; % second derivative of the refractive index with respect to omega
refr_omega_ppp=diff(refr_omega_pp)/delta_omega; % third derivative of the refractive index with respect to omega

refr_omega0=sellmeier_FS(c./(omega0/2/pi).*1000); % refractive index at omega0;
refr_omega0_p=interp1(omega_p,refr_omega_p,omega0); % first derivative of the refractive index at omega0;
refr_omega0_pp=interp1(omega_pp,refr_omega_pp,omega0); % ...
refr_omega0_ppp=interp1(omega_ppp,refr_omega_ppp,omega0);

apex_rad=apex_angle/180.*pi;
theta1=asin(sin(theta0)./refr_omega0);
theta3=asin(refr_omega0.*sin(apex_rad-asin(sin(theta0)./refr_omega0)));
theta1_p=-tan(theta1)./refr_omega0.*refr_omega0_p; % first derivative of theta1 with respect to omega
theta1_pp=tan(theta1)./refr_omega0^2.*(1+1./cos(theta1)^2).*refr_omega0_p^2-tan(theta1)./refr_omega0.*refr_omega0_pp; %  second derivative of theta1 with respect to omega
theta3_p=sin(apex_rad)./cos(theta1)./cos(theta3).*refr_omega0_p; % first derivative of theta3 with respect to omega
theta3_pp=sin(apex_rad)./cos(theta1)./cos(theta3).*((sin(apex_rad).*tan(theta3)./cos(theta1)./cos(theta3)-tan(theta1)^2./refr_omega0).*refr_omega0_p^2+refr_omega0_pp);
Lg=(t-s.*tan(theta3)).*sin(apex_rad)./cos(theta1);
Lg_p=Lg.*tan(theta1).*theta1_p-s.*sin(apex_rad)./cos(theta1)./cos(theta3)^2.*theta3_p; % first derivative of Lg with respect to omega

% PsiN = (d^N Psi)/(d omega^N), where Psi is the phase as a function of omega

Psi0=omega0/c*(refr_omega0.*Lg+s./cos(theta3)-sin(theta0).*((t-s.*tan(theta3)).*(cos(apex_rad)+sin(apex_rad).*tan(theta1))-a));
Psi1=Psi0./omega0+omega0/c.*Lg.*refr_omega0_p;
Psi2=Lg./c.*(2.*refr_omega0_p+omega0.*refr_omega0_pp)-omega0./c.*(refr_omega0.*Lg.*theta1_p^2+s./cos(theta3).*theta3_p^2); % equation 2.101 on p114
Psi3_first=1./c.*Lg_p.*(2.*refr_omega0_p+omega0.*refr_omega0_pp)+Lg./c.*(3.*refr_omega0_pp+omega0.*refr_omega0_ppp);
Psi3_second=-1./c.*(refr_omega0.*Lg.*theta1_p^2+s./cos(theta3).*theta3_p^2)-omega0./c.*(refr_omega0_p.*Lg.*theta1_p^2+refr_omega0.*Lg_p.*theta1_p^2+2.*refr_omega0.*Lg.*theta1_p.*theta1_pp+s.*tan(theta3)./cos(theta3).*theta3_p^3+2.*s./cos(theta3).*theta3_pp.*theta3_p);
Psi3=Psi3_first+Psi3_second;

varargout{1}=Psi0;
varargout{2}=Psi1;
varargout{3}=Psi2;
varargout{4}=Psi3;
end
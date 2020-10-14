clear;
center_wl=0.78; % [micron]
% apex_dist=800; % OO' distance in [mm]
travel_dist=155+525; % BB' distance in [mm]
apex_angle=67.88;
apex_rad=apex_angle/180*pi;
n0=sellmeier_FS(center_wl);
brewster=atan(n0);
mindev_angle=asin(n0*sin(apex_rad/2));
theta0=55.5; %[degrees]
theta0_rad=theta0/180*pi;
theta1=asin(sin(theta0_rad)/n0);
theta3=asin(n0.*sin(apex_rad-asin(sin(theta0_rad)./n0)));
omega=2*pi*(300:10:750)/1000;
omega0=300./(center_wl*1000)*2*pi;
% apex_angle=pi-2*theta0;
a=3; % distance from 1st prism apex to where the beam hits the first interface
insertion1=a*cos(apex_rad/2); % 1st prism insertion in [mm], measured along 'y', the coordinate along the angle bisector of the apex angle
insertion2=4; % 2nd prism insertion in [mm], measured along 'y'.
b=a*cos(apex_rad)+a*sin(apex_rad)*tan(apex_rad-theta1); % distance from 1st prism apex to where the beam hits the second interface (= OB distance)
t0=travel_dist*sin(theta3)+b;
s0=travel_dist*cos(theta3);
M=[cos(apex_rad/2) -sin(apex_rad/2); sin(apex_rad/2) cos(apex_rad/2)];
Minv=inv(M);
for ind1=1:length(travel_dist)
    temp=M*[t0(ind1); s0(ind1)];
    y0(ind1)=temp(1);
    x0(ind1)=temp(2);
    for ind2=1:length(insertion2)
        y(ind1+ind2-1)=y0(ind1)+insertion2(ind2);
        temp=Minv*[y(ind1+ind2-1); x0(ind1)];
        t(ind1+ind2-1)=temp(1);
        s(ind1+ind2-1)=temp(2);
    end
end

[phi0 phi1 phi2 phi3]=prism_compressor(center_wl,apex_angle,theta0,a,t,s);
% end
%%
str_title{1}=['\lambda_0= ' num2str(center_wl*1000) ' nm; apex angle = ' num2str(apex_angle) '^\circ; FS prisms'];
str_title{3}=['\theta_0: ' num2str(roundP(theta0/pi*180,2)) '^\circ; \theta_3: ' num2str(roundP(theta3/pi*180,2)) '^\circ; minimum dev. angle: ' num2str(roundP(mindev_angle/pi*180,2)) '^\circ'];
str_title{2}=['1st prism insertion = ' num2str(insertion1/10) ' cm; 2nd prism insertion = ' num2str(insertion2/10) ' cm']
figure;
subplot(211);plot(travel_dist,2*phi2,'k-') % factor of 2 is inserted for double pass geometry
hold on;plot(travel_dist,zeros(size(travel_dist)),'r--')
ylabel('GDD [fs^2]')
setfigP
title(str_title)
subplot(212);plot(travel_dist,2*phi3,'k-')
hold on;plot(travel_dist,zeros(size(travel_dist)),'r--')
ylabel('TOD [fs^3]')
xlabel('prism distance (BB'') [mm]')
setfigP

for ind1=1:length(omega0)
    [phi2(ind1) phi3(ind1)]=prism_compressor(omega,omega0(ind1),apex_angle,theta0,t,s);
end
str_title{1}=['insertion= ' num2str(insertion) ' mm; apex angle = ' num2str(apex_angle) ' [deg]; apex distance = ' num2str(apex_dist/10) ' cm'];
str_title{2}=['\theta_0= ' num2str(roundP(theta0/pi*180,2)) ' [deg]'];
% 
% figure;
% subplot(211);plot(center_wl*1000,phi2)
% ylabel('GDD [fs^2]')
% title(str_title)
% subplot(212);plot(center_wl*1000,phi3)
% ylabel('TOD [fs^3]')
% xlabel('\lambda_0 [nm]')
%%
figure;plotyyP2(insertion2,map2colvec(phi2),map2colvec(phi3))
title(str_title)
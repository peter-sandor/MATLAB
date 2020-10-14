% Framework, notation and equations are from the book
% Ultrashort Laser Pulse Phenomena by Diels and Rudolph

clear;
lambda0=0.795; % reference wavelength in [micron]
% apex_dist=800; % OO' distance in [mm]
travel_dist_vec=100:20:1000; % BB' distance in [mm]
N=length(travel_dist_vec);
% apex_angle=56.3504; % apex angle for the two prisms in [degrees]
prism_material='SF10';
apex_angle=60;
apex_rad=apex_angle/180*pi;

[center_wl2,n0]=sellmeier(lambda0,prism_material);
theta_mindev_rad=asin(n0*sin(apex_rad/2)); % incident angle for minimum deviation
theta_mindev=theta_mindev_rad/pi*180;
theta0=theta_mindev; % incident angle is the one for minimum deviation
% theta0=65; % incident angle in [degrees]
theta0_rad=theta0/180*pi;
phi0=zeros([N 1]);
phi1=zeros([N 1]);
phi2=zeros([N 1]);
phi3=zeros([N 1]);

for ind3=1:N
    travel_dist=travel_dist_vec(ind3);
    center_wl=lambda0;%(ind3);
    theta1_rad=asin(sin(theta0_rad)/n0);
    theta3_rad=asin(n0.*sin(apex_rad-asin(sin(theta0)./n0)));
%     omega=2*pi*(300:10:750)/1000;
%     omega0=300./(center_wl*1000)*2*pi;
%     apex_angle=pi-2*theta0;
    a=43; % distance from 1st prism apex to where the beam hits the first interface (OA)
    b2=50; % distance from 2nd prism apex to where the beam hits the first interface (O'B')
    insertion1=a*cos(apex_rad/2); % 1st prism insertion in [mm], measured along 'y', the coordinate along the angle bisector of the apex angle
    insertion2=b2*cos(apex_rad/2); % 2nd prism insertion in [mm], measured along 'y'.
    b=a*cos(apex_rad)+a*sin(apex_rad)*tan(apex_rad-theta1_rad); % distance from 1st prism apex to where the beam hits the second interface (= OB distance)
    a2=b2*cos(apex_rad)+b2*sin(apex_rad)*tan(theta1_rad);
    t0=travel_dist*sin(theta3_rad)+b+b2;
    s0=travel_dist*cos(theta3_rad);
    % Matrix M is for the transformation between the rectangular
    % coordinate systems [x y] and [s t].
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

    [phi0(ind3) phi1(ind3) phi2(ind3) phi3(ind3)]=prism_compressor(center_wl,apex_angle,theta0,a,t,s);
end
%%
str_title{1}=['\lambda_0= ' num2str(center_wl*1000) ' nm; apex angle = ' num2str(apex_angle) '^\circ; ' prism_material ' prisms'];
str_title{3}=['\theta_0: ' num2str(roundP(theta0,2)) '^\circ; \theta_3: ' num2str(roundP(theta3_rad/pi*180,2)) '^\circ; minimum dev. angle: ' num2str(roundP(theta_mindev,2)) '^\circ'];
str_title{2}=['1st prism insertion = ' num2str(insertion1) ' mm; 2nd prism insertion = ' num2str(insertion2) ' mm'];
figure;
subplot(211);plot(travel_dist_vec,2*phi2,'k-') % factor of 2 is inserted for double pass geometry
hold on;plot(travel_dist,zeros(size(travel_dist)),'r--')
ylabel('GDD [fs^2]')
setfigP
title(str_title)
subplot(212);plot(travel_dist_vec,2*phi3,'k-')
hold on;plot(travel_dist,zeros(size(travel_dist)),'r--')
ylabel('TOD [fs^3]')
% xlabel('wavelength [micron]')
xlabel('prism distance (BB'') [mm]')
setfigP

% for ind1=1:length(omega0)
%     [phi2(ind1) phi3(ind1)]=prism_compressor(omega,omega0(ind1),apex_angle,theta0,t,s);
% end
% str_title{1}=['insertion= ' num2str(insertion) ' mm; apex angle = ' num2str(apex_angle) ' [deg]; apex distance = ' num2str(apex_dist/10) ' cm'];
% str_title{2}=['\theta_0= ' num2str(roundP(theta0/pi*180,2)) ' [deg]'];
% 
% figure;
% subplot(211);plot(center_wl*1000,phi2)
% ylabel('GDD [fs^2]')
% title(str_title)
% subplot(212);plot(center_wl*1000,phi3)
% ylabel('TOD [fs^3]')
% xlabel('\lambda_0 [nm]')
%%
N=2;
xprism=y0*tan(apex_rad/2);
framex=map2colvec(-xprism:(x0+2*xprism)/(N-1):x0+xprism);
framey=map2colvec(0:y0/(N-1):y0);
line_prism1a=[map2colvec(-xprism:tan(apex_rad/2)*y0/(N-1):0), map2colvec(y0:-y0/(N-1):0)];
line_prism1b=[map2colvec(0:tan(apex_rad/2)*y0/(N-1):tan(apex_rad/2)*y0), map2colvec(0:y0/(N-1):y0)];
line_prism2a=[map2colvec(x0-xprism:tan(apex_rad/2)*y0/(N-1):x0), map2colvec(0:y0/(N-1):y0)];
line_prism2b=[map2colvec(x0:tan(apex_rad/2)*y0/(N-1):x0+xprism), map2colvec(y0:-y0/(N-1):0)];
line_beam1=[map2colvec(-xprism:(xprism-a*sin(apex_rad/2))/(N-1):-a*sin(apex_rad/2)), map2colvec(a*cos(apex_rad/2)+(xprism-a*sin(apex_rad/2))*tan(theta0_rad-apex_rad/2):-(xprism-a*sin(apex_rad/2))*tan(theta0_rad-apex_rad/2)/(N-1):a*cos(apex_rad/2))];
line_beam2=[map2colvec(-a*sin(apex_rad/2):(a+b)*sin(apex_rad/2)/(N-1):b*sin(apex_rad/2)), map2colvec(a*cos(apex_rad/2):(b-a)*cos(apex_rad/2)/(N-1):b*cos(apex_rad/2))];
line_beam3=[map2colvec(b*sin(apex_rad/2):travel_dist*cos(theta3_rad-apex_rad/2)/(N-1):b*sin(apex_rad/2)+travel_dist*cos(theta3_rad-apex_rad/2)), map2colvec(b*cos(apex_rad/2):travel_dist*sin(theta3_rad-apex_rad/2)/(N-1):b*cos(apex_rad/2)+travel_dist*sin(theta3_rad-apex_rad/2))];
line_beam4=[map2colvec(b*sin(apex_rad/2)+travel_dist*cos(theta3_rad-apex_rad/2):(a2+b2)*sin(apex_rad/2)/(N-1):b*sin(apex_rad/2)+travel_dist*cos(theta3_rad-apex_rad/2)+(a2+b2)*sin(apex_rad/2)), map2colvec(b*cos(apex_rad/2)+travel_dist*sin(theta3_rad-apex_rad/2):(b2-a2)*cos(apex_rad/2)/(N-1):b*cos(apex_rad/2)+travel_dist*sin(theta3_rad-apex_rad/2)+(b2-a2)*cos(apex_rad/2))];
line_beam5=[map2colvec(b*sin(apex_rad/2)+travel_dist*cos(theta3_rad-apex_rad/2)+(a2+b2)*sin(apex_rad/2):(max(framex)-(b*sin(apex_rad/2)+travel_dist*cos(theta3_rad-apex_rad/2)+(a2+b2)*sin(apex_rad/2)))/(N-1):max(framex)), map2colvec(b*cos(apex_rad/2)+travel_dist*sin(theta3_rad-apex_rad/2)+(b2-a2)*cos(apex_rad/2):-(a2*cos(apex_rad/2)+(xprism-a2*sin(apex_rad/2))*tan(theta0_rad-apex_rad/2))/(N-1):b*cos(apex_rad/2)+travel_dist*sin(theta3_rad-apex_rad/2)+(b2-a2)*cos(apex_rad/2) - (a2*cos(apex_rad/2)+(xprism-a2*sin(apex_rad/2))*tan(theta0_rad-apex_rad/2)))];

dfigure;axes;
line(line_prism1a(:,1),line_prism1a(:,2),'color','k')
line(line_prism1b(:,1),line_prism1b(:,2),'color','k')
line(line_prism2a(:,1),line_prism2a(:,2),'color','k')
line(line_prism2b(:,1),line_prism2b(:,2),'color','k')

line(line_beam1(:,1),line_beam1(:,2),'color','b')
line(line_beam2(:,1),line_beam2(:,2),'color','b')
line(line_beam3(:,1),line_beam3(:,2),'color','b')
line(line_beam4(:,1),line_beam4(:,2),'color','b')
line(line_beam5(:,1),line_beam5(:,2),'color','b')

line(framex,zeros(size(framex)),'color','r','linestyle','--')
line(framex,y0*ones(size(framex)),'color','r','linestyle','--')
line(min(framex)*ones(size(framey)),framey,'color','r','linestyle','--')
line(max(framex)*ones(size(framey)),framey,'color','r','linestyle','--')
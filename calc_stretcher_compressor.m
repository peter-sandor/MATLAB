function [ deltaphi2,deltaphi3 ] = calc_stretcher_compressor(b,theta,L,gamma,lambda)
% mm and fs are used
% b= perpendicular compressor grating separation in mm
% theta= compressor incident angle in degrees
% L= stretcher l_eff in mm
% gamma= stretcher incident angle in degrees

N=[length(b) length(theta) length(L) length(gamma)];

c=3*10^(-4);
n=1200; % lines/mm
% lambda=785*10^-6;
m=1;
phi2=zeros(N);
phi3=zeros(N);
phi2s=zeros(N);
phi3s=zeros(N);
gamma1=zeros([N(4) 1]);
B=zeros([N(3) N(4)]);

for ind1=1:N(1)
    for ind2=1:N(2)
        for ind3=1:N(3)
            for ind4=1:N(4)
                phi2(ind1,ind2,ind3,ind4)=-b(ind1).*lambda./(pi.*c^2).*(lambda.*n)^2./(1-(sin(theta(ind2).*pi./180)-lambda.*n).^2).^(3/2);
                phi3(ind1,ind2,ind3,ind4)=-3.*lambda./(2.*pi.*c).*phi2(ind1,ind2,ind3,ind4).*(1+lambda.*n.*(lambda.*n-sin(theta(ind2).*pi./180))./(1-(lambda.*n-sin(theta(ind2).*pi./180)).^2));
%                 p2=-b(ind1)*n^2*lambda^3/(pi*c^2*cos(18.63*pi/180)^3);
                gamma1(ind4)=abs(asin(sin(gamma(ind4).*pi/180)-lambda.*n).*180/pi);  % use grating equation to calculate diffracted angle from incident angle
                B(ind3,ind4)=-L(ind3)*2*cos(gamma1(ind4)*pi/180); % calculate perpendicular grating separation for the stretcher
                phi2s(ind1,ind2,ind3,ind4)=-B(ind3,ind4).*lambda./(pi.*c^2).*(lambda.*n)^2./(1-(sin(gamma(ind4).*pi./180)-lambda.*n).^2).^(3/2);
                phi3s(ind1,ind2,ind3,ind4)=-3.*lambda./(2.*pi.*c).*phi2s(ind1,ind2,ind3,ind4).*(1+lambda.*n.*(lambda.*n-sin(gamma(ind4).*pi./180))./(1-(lambda.*n-sin(gamma(ind4).*pi./180)).^2));
            end
        end
    end
end
deltaphi2= squeeze(phi2+phi2s); % net second order phase in [fs^2]
deltaphi3= squeeze(phi3+phi3s); % net third order phase in [fs^3]
end
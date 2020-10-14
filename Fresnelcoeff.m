function [Rp,Rs]=Fresnelcoeff(theta,n1,n2)
% input: theta in [deg]
theta=theta/180*pi;
Rp=abs((n1*sqrt(1-(n1/n2*sin(theta)).^2)-n2*cos(theta))./(n1*sqrt(1-(n1/n2*sin(theta)).^2)+n2*cos(theta))).^2;
Rs=abs((n1*cos(theta)-(n2*sqrt(1-(n1/n2*sin(theta)).^2)))./(n2*sqrt(1-(n1/n2*sin(theta)).^2)+n1*cos(theta))).^2;
end
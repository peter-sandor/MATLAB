l_side = 4;
lambda = 0.3:0.01:1;
n = sellmeier(lambda,'Sapphire_o');
figure;plot(lambda,n,'k')
n_Sapph = sellmeier(0.8,'Sapphire_o');
% thetaB = asin(sqrt(n_Sapph^2/(1+n_Sapph^2)))/pi*180;
thetaB = atan(n_Sapph)/pi*180; % [deg]
beta = asin(1/n_Sapph*sin(thetaB/180*pi))/pi*180; % [deg]
l_path = (l_side*sin(thetaB/180*pi))/cos(beta/180*pi);

-log(0.08)/0.4/n_Sapph
exp(-6.3*0.4)
exp(-4.5*0.4*n_Sapph)
function [rho,f0] = invAbel5(x,F0,N)

% 1D Abel inversion algorithm based on direct evaluation of the inversion
% integral (see [http://mathworld.wolfram.com/AbelTransform.html] and [http://en.wikipedia.org/wiki/Abel_transform])
%
% Input Arguments:
% x     --  abscissa vector for the input vector F0 (x and F0 has to be the same length)
% F0    --  function to be inverted
% N     --  desired length for the output vectors rho and f0
%
% Output Arguments:
% rho   --  abscissa vector for the output vector f0
% f0    --  result of the Abel-inversion

Nx=length(x);
N_samp=10*Nx;
x_min=min(x);
x_max=max(x);
x_step=(x_max-x_min)/(Nx);

eps2=1e-6;
param1=3.5; % Sets offset between abscissa vectors x and rho; Magic number: 3.5
scale1=1;
eps=x_step/20;
func = @(x) scale1*(x).^1; % this one could be used to generate abscissa vectors with nonequidistant spacing

x2=func((x_min:(x_max-x_min)/(N_samp-1):x_max)')/(func(x_max))*x_max;
rho1=x2;
x2=x2+(x2(2)-x2(1))/param1;
eps=(x2(2)-x2(1))/param1/20;
x2=x2(abs(x2)<=x_max);

F1=interp1(x,F0,x2,'linear');
rho1=rho1(rho1<max(x2));
Nrho1=length(rho1);

rho_min=min(rho1);
rho_max=max(rho1);
rho_step=(rho_max-rho_min)/(N-1);
rho=(rho_min:rho_step:rho_max)';

dF2=interp1(x,gradient(F0,x),x2);
dF=gradient(F1,x2);

% func_invtd = zeros([Nrho1,1]);
func_invtd2 = zeros([Nrho1,1]);
for k=1:Nrho1
%     func_invtd(k)=-1/pi*trapz(x2((x2-rho1(k))>=eps),dF((x2-rho1(k))>=eps).*1./sqrt(x2((x2-rho1(k))>=eps).^2-(rho1(k))^2));
    func_invtd2(k)=-1/pi*sum(dF((x2-rho1(k))>=eps).*1./sqrt(x2((x2-rho1(k))>=eps).^2-(rho1(k))^2).*gradient(x2((x2-rho1(k))>=eps)));
end
f0 = interp1(rho1,func_invtd2,rho,'linear');
end
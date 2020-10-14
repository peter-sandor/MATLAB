function [x,F0] = Abel(rho,f0,N)

Nrho=length(rho);
factor=10;
N_samp=factor*Nrho;
rho_min=min(rho);
rho_max=max(rho);
rho_step=(rho_max-rho_min)/(Nrho);

scale1=1;
param1=-0.7; % magic value: -0.7
func = @(x) scale1*(x).^1; % this one could be used to generate abscissa vectors with nonequidistant spacing
rho2=func((rho_min:(rho_max-rho_min)/(N_samp-1):rho_max)')/(func(rho_max))*rho_max;
x1=rho2;
rho2=rho2+(rho2(2)-rho2(1))*param1;
eps=(rho2(2)-rho2(1))/param1/20;
rho2=rho2(abs(rho2)<=rho_max);

f1=interp1(rho,f0,rho2,'linear');
x1=x1(x1<max(rho2));
Nx=length(x1);

x_min=min(x1);
x_max=max(x1);
x_step=(x_max-x_min)/(N-1);
x=(x_min:x_step:x_max)';

func_projd = zeros([Nx 1]);
for k=1:Nx
    func_projd(k)=2*sum(f1(x1>x1(k)).*rho2(x1>x1(k))./sqrt(rho2(x1>x1(k)).^2-x1(k).^2)*rho_step/factor);
end
F0 = interp1(x1,func_projd,x,'linear');

end
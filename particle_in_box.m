h = 6.636e-34; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
L = [1e-8 1e-9 1e-10]; % [m]
eV = 1.602e-19; % 1 J in electronvolts
n = 1;
E = h^2./(8*L.^2*mass)*n.^2/eV; % [eV]
figure;plot(L,E,'k.')

%%
n = (1:7)';
a = 1/(sqrt(2)*pi).*(1./(n./2-1).*sin((n./2-1)*pi)-1./(n./2+1).*sin((n./2+1)*pi));

a(2) = 1/sqrt(2);
for k = 4:2:20
a(k)= 0;
end
%%
h = 6.636e-34; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
L = 1e-10; % [m]
x = 0:0.01:1;
amp = zeros([length(x),length(n)]);
omega = h/(2*pi)*pi^2*n.^2/(8*mass*L);
T = 2*pi/omega(1);
t = T/2;
for k = 1:length(n)
amp(:,k) = a(k)*sin(n(k)*pi*x)*real(exp(-i*omega(k)*t));
end

psi = zeros([length(x),1]);
psi = sum(amp(:,1:7),2);

figure;plot(x,amp(:,1),'y')
hold on; plot(x,amp(:,2), 'r')
plot(x,amp(:,3), 'b')
plot(x,amp(:,5), 'g')
plot(x,amp(:,7), 'm')
plot(x,psi,'k.')
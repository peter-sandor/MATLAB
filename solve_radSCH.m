%% Solve nonrelativistic radial Schroedinger equation in atomic units

zmax = 50; % units of length are Bohr = 0.53 Angstrom
zmin = 0;
N = 100;
zstp = (zmax-zmin)/(N-1);
z = (zmin:zstp:zmax);
n=20; % principle quantum number; plugged into the expression for energy, in Rydberg (1 Ry = 27.2 eV)
k=1; % k=Z/(4*pi*eps0)=Z in atomic units
l=0; % orbital angular momentum quantum number

% initial conditions [x(1);x(2);x(3)]
init_cond = [zmax,0,1];
func_radSCH = @(z,x) [1; x(3); -(-1/n^2-2*k/x(1)+l*(l+1)/x(1)^2)*x(2)];
[zout, yout] = ode45(func_radSCH,[zmax zmin],init_cond);
R=abs(yout(1:end-2,2)./zout(1:end-2)).^1/sum(abs(yout(1:end-2,2)).^2,1);
dfigure;
plot(zout,yout(:,2).^2,'k')
% hold on; plot(zout(1:end-2),R,'r')
xlim([0 zmax])
xlabel('r [Bohr]')
ylabel(['\chi_{' num2str(n) ',' num2str(l) '}(r)'])
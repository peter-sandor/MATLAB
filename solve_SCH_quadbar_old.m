syms x t A alphat xt pt gammat smi hbar mass omega diffxt diffpt diffalphat diffgammat;

psi = A*exp(-alphat*(x-xt)^2+smi/hbar*pt*(x-xt)+smi/hbar*gammat);

lhs = smi*hbar*(diff(psi,xt)*diffxt + diff(psi,pt)*diffpt + diff(psi,alphat)*diffalphat + diff(psi,gammat)*diffgammat);
rhs = -hbar^2/(2*mass)*(diff(psi,'x',2)-1/2*mass*omega^2*x^2);

pretty(lhs)
pretty(rhs)
%%
omega = 2*pi*7.5e8;
alpha_0 = 1e12; 
% time_scale=2*pi/omega; % harmonic oscillation scenario
% time_scale = 4e-10; % backward scattering scenario
time_scale = 5e-11; % crossing the barrier scenario
N_delay = 64;
delay_rng = 2*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
% vec_0 = [3e-7;0;1.0395e+012;1e-34;]; % "free fall"
% vec_0 = [-1e-4;4e-25;alpha_0;1e-34;]; % "backward scattering"
vec_0 = [-1e-4;1e-24;alpha_0;1e-34;]; % "backward scattering"

[T, SCH_qb_ode45] = ode45(@func_SCH_quadbar, delay, vec_0);

% figure;
% subplot(211)
% plot(T,SCH_qb_ode45(:,1)/max(SCH_qb_ode45(:,1)),'k')
% hold on
% plot(T,SCH_qb_ode45(:,2)/max(SCH_qb_ode45(:,2)),'r')
% subplot(212)
% plot(T,abs(SCH_qb_ode45(:,3))/max(abs(SCH_qb_ode45(:,3))),'b')
% hold on
% plot(T,abs(SCH_qb_ode45(:,4))/max(abs(SCH_qb_ode45(:,4))),'g')

%%
Amp = ((2*alpha_0)/pi)^(1/4);
h_bar = 6.626e-34;
N_x = 256;
% N_padded = 2048;

% Physical distance range for the "backward scattering" scenario
% x_rng = 1e-4;
% x_stp = x_rng/N_x;
% x = (-x_rng:x_stp:0)'; % [m]

% Physical distance range for the "crossing the barrier" scenario
x_rng = 1e-5;
x_stp = x_rng/N_x;
x = ((-9*x_rng/10+x_stp):x_stp:x_rng/10)'; % [m]

% Physical distance range for the "free fall" scenario
% x_rng = 1e-4;
% x_stp = x_rng/N_x;
% x = (0:x_stp:x_rng)'; % [m]

h_fig = figure;
for k = 1:length(delay)
%     Amp = (2*abs(SCH_qb_ode45(k,3))/pi)^(1/4);
    Psi = Amp.*exp(-SCH_qb_ode45(k,3).*(x-SCH_qb_ode45(k,1)).^2+i/h_bar*SCH_qb_ode45(k,2)*(x-SCH_qb_ode45(k,1))+i/h_bar.*SCH_qb_ode45(k,4));
    plot(x,abs(Psi),'k.')
    % ylim([-4 4])
    xlabel('Physical distance [m]')
%     axis([min(x),max(x),-1e4,5e4])
    hold on;
    plot(x,imag(Psi),'b')
    plot(x,real(Psi),'r')
    hold off;
    F(k) = getframe(h_fig);
end
close(h_fig)
filename = 'gaussian_at_quadratic_barrier_4.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
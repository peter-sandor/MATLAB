syms x t A alphat xt pt gammat smi hbar mass omega diffxt diffpt diffalphat diffgammat;

psi = A*exp(-alphat*(x-xt)^2+smi/hbar*pt*(x-xt)+smi/hbar*gammat);

lhs = smi*hbar*(diff(psi,xt)*diffxt + diff(psi,pt)*diffpt + diff(psi,alphat)*diffalphat + diff(psi,gammat)*diffgammat);
rhs = -hbar^2/(2*mass)*(diff(psi,'x',2)-1/2*mass*omega^2*x^2);

pretty(lhs)
pretty(rhs)
%%
omega = 2*pi*7.5e8;
time_scale=2*pi/omega;
N_delay = 64;
delay_rng = 2*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
vec_0 = [3e-7;0;1.0395e+012;1e-34;];
% vec_0 = [1;0.1;0;0];
% 
[T, SCH_qp_ode45] = ode45(@func_SCH_quadpot, delay, vec_0);

figure;
subplot(211)
plot(T,SCH_qp_ode45(:,1)/max(SCH_qp_ode45(:,1)),'k')
hold on
plot(T,SCH_qp_ode45(:,2)/max(SCH_qp_ode45(:,2)),'r')
subplot(212)
plot(T,abs(SCH_qp_ode45(:,3))/max(abs(SCH_qp_ode45(:,3))),'b')
hold on
plot(T,abs(SCH_qp_ode45(:,4))/max(abs(SCH_qp_ode45(:,4))),'g')

%%
Amp = 1;
h_bar = 6.626e-34;
N_x = 256;
% N_padded = 2048;
x_rng = 5e-6;
x_stp = x_rng/N_x;
x = ((-x_rng/2+x_stp):x_stp:x_rng/2)'; % [m]
h_fig = figure;
for k = 1:length(delay)
    Psi = Amp.*exp(-SCH_qp_ode45(k,3).*(x-SCH_qp_ode45(k,1)).^2+i/h_bar*SCH_qp_ode45(k,2)*(x-SCH_qp_ode45(k,1))+i/h_bar.*SCH_qp_ode45(k,4));
    plot(x,abs(Psi),'k.')
    xlabel('Physical distance [m]')
%     axis([min(x),max(x),-1e4,5e4])
    hold on;
    plot(x,imag(Psi),'b')
    plot(x,real(Psi),'r')
    hold off;
    F(k) = getframe(h_fig);
end
close(h_fig)
filename = 'test.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
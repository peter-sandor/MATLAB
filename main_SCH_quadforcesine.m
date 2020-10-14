
omega = 2*pi*7.5e8;
time_scale=2*pi/omega;
N_delay = 64;
delay_rng = 4*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
% initial conditions for the 4 coupled variables, which are:
% [x_t; p_t; alpha_t, gamma_t]
init_cond = [3e-7;0;1.0395e+012;1e-34;];
% parameters: [electron mass [kg]; Planck-constant/(2*pi) [Js]; Natural frequency [Hz]; Natural frequency of driving force [Hz]; Phase offset of driving force; Driving force amplitude [N]]
param = [9e-31;6.626e-34/2/pi;2*pi*7.5e8; 0 *2*pi*7.5e8;pi/2;1e-17];

[T,SCH_qp_ode45] = solve_SCH_quadforcesine(delay,init_cond,param);

figure;
subplot(211)
plot(T,SCH_qp_ode45(:,1)/max(abs(SCH_qp_ode45(:,1))),'k')
hold on
plot(T,SCH_qp_ode45(:,2)/max(abs(SCH_qp_ode45(:,2))),'r')
subplot(212)
plot(T,abs(SCH_qp_ode45(:,3))/max(abs(SCH_qp_ode45(:,3))),'b')
hold on
plot(T,abs(SCH_qp_ode45(:,4))/max(abs(SCH_qp_ode45(:,4))),'g')

%%
Amp = 1;
h_bar = 6.626e-34;
N_x = 256;
% N_padded = 2048;
x_rng = 8*max(SCH_qp_ode45(:,1))-min(SCH_qp_ode45(:,1));
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
filename = 'teszt.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
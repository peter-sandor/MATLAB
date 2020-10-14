mass = 9e-31; % [kg]
omega = 2*pi*7.5e8;
h_bar = 6.626e-34; % [Js]
time_scale=2*pi/omega;
param_a = mass*omega/(2*h_bar);
alpha_0 = param_a;
N_delay = 64;
delay_rng = 2*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
% initial conditions for the 4 coupled variables, which are:
% [x_t; p_t; alpha_t, gamma_t]
init_cond = [-3e-7;0;alpha_0;1e-34;];
% parameters: [electron mass [kg]; Planck-constant/(2*pi) [Js]; Natural frequency [Hz]; Natural frequency of driving force [Hz]; Phase offset of driving force; Driving force amplitude [N]]
param = [mass;h_bar/2/pi;omega;];

[T, SCH_qp_ode45] = solve_SCH_quadpot(delay,init_cond,param);

figure;
subplot(311)
plot(delay,SCH_qp_ode45(:,1),'k')
subplot(312)
plot(delay,SCH_qp_ode45(:,2),'r')
subplot(313)
plot(delay,abs(SCH_qp_ode45(:,3)),'b')

% figure;
% subplot(211)
% plot(T,SCH_qp_ode45(:,1)/max(SCH_qp_ode45(:,1)),'k')
% hold on
% plot(T,SCH_qp_ode45(:,2)/max(SCH_qp_ode45(:,2)),'r')
% subplot(212)
% plot(T,abs(SCH_qp_ode45(:,3))/max(abs(SCH_qp_ode45(:,3))),'b')
% hold on
% plot(T,abs(SCH_qp_ode45(:,4))/max(abs(SCH_qp_ode45(:,4))),'g')

%%
Amp = 1;
N_x = 256;
% N_padded = 2048;
x_rng = 2e-6;
x_stp = x_rng/N_x;
x = ((-x_rng/2+x_stp):x_stp:x_rng/2); % [m]
p_rng = 5e-27;
p_stp = p_rng/N_x;
p = ((-p_rng/2+p_stp):p_stp:p_rng/2)';
h_fig = figure;
for k = 1:length(delay)
    PSI_WR1 = exp(-(p-SCH_qp_ode45(k,2)).^2/(h_bar^2*2*real(SCH_qp_ode45(k,3))));
    PSI_WR2 = exp(-2*abs(SCH_qp_ode45(k,3))^2/(real(SCH_qp_ode45(k,3)))*(x-SCH_qp_ode45(k,1)).^2);
    PSI_WR3 = exp(-2*imag(SCH_qp_ode45(k,3))/(h_bar*real(SCH_qp_ode45(k,3)))*(p-SCH_qp_ode45(k,2))*(x-SCH_qp_ode45(k,1)));
    Psi_WR = Amp^2/(sqrt(2*pi*real(SCH_qp_ode45(k,3)))*h_bar*2*pi).*PSI_WR1*PSI_WR2.*PSI_WR3;
%     Psi_WR = Amp^2/(sqrt(2*pi*real(alpha_t(k)))*h_bar*2*pi)*exp(-2*abs(alpha_t(k))^2/(real(alpha_t(k))).*(x'-x_t(k)).^2)*exp(-(p'-p_t(k)).^2/(h_bar^2*2*real(alpha_t(k))))*exp(+2*imag(alpha_t(k))/(h_bar*real(alpha_t(k)))*(x'-x_t(k))*(p'-p_t(k)));
    imagesc(x,p,Psi_WR)
    
    xlabel('Physical distance [m]')
    ylabel('Momentum [kg*m/s]')
    F(k) = getframe(h_fig);
end
close(h_fig)
filename = 'test.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
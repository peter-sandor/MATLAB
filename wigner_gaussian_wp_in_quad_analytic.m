clear;
h_bar = 6.636e-34/2/pi; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
omega = 2*pi*7.5e8;
param_a = mass*omega/(2*h_bar);
% alpha_0 = 1e12; 
alpha_0 = param_a*5;

% alpha_0 = param_a;
x_0 = -3e-7;
% p_0 = mass*omega*5e-8;
p_0 = 0;
Amp = (2*alpha_0/pi)^(1/4);
time_scale=2*pi/omega; % [s]
N_x = 256;
% N_padded = 2048;
x_rng = 2e-6;
x_stp = x_rng/N_x;
x = ((-x_rng/2+x_stp):x_stp:x_rng/2); % [m]
p_rng = 5e-27;
p_stp = p_rng/N_x;
p = ((-p_rng/2+p_stp):p_stp:p_rng/2)';
% p = flipdim(((-p_rng/2+p_stp):p_stp:p_rng/2)',1); % [kg*m/s]

V_quad = 1/2*mass*omega^2*x.^2;
N_delay = 64;
delay_rng = 2*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
x_t = x_0*cos(omega*delay)+p_0/(mass*omega)*sin(omega*delay);
p_t = p_0*cos(omega*delay)-mass*omega*x_0*sin(omega*delay);
alpha_t = param_a*(alpha_0*cos(omega*delay)+i*param_a*sin(omega*delay))./(param_a*cos(omega*delay)+i*alpha_0*sin(omega*delay));

h_fig = figure;
for k = 1:length(delay)
% k = 64;
    t=delay(k); % [s]
    PSI_WR1 = exp(-(p-p_t(k)).^2/(h_bar^2*2*real(alpha_t(k))));
    PSI_WR2 = exp(-2*abs(alpha_t(k))^2/(real(alpha_t(k)))*(x-x_t(k)).^2);
    PSI_WR3 = exp(-2*imag(alpha_t(k))/(h_bar*real(alpha_t(k)))*(p-p_t(k))*(x-x_t(k)));
    Psi_WR = Amp^2/(sqrt(2*pi*real(alpha_t(k)))*h_bar*2*pi).*PSI_WR1*PSI_WR2.*PSI_WR3;
%     Psi_WR = Amp^2/(sqrt(2*pi*real(alpha_t(k)))*h_bar*2*pi)*exp(-2*abs(alpha_t(k))^2/(real(alpha_t(k))).*(x'-x_t(k)).^2)*exp(-(p'-p_t(k)).^2/(h_bar^2*2*real(alpha_t(k))))*exp(+2*imag(alpha_t(k))/(h_bar*real(alpha_t(k)))*(x'-x_t(k))*(p'-p_t(k)));
    imagesc(x,p,Psi_WR)
    
    xlabel('Physical distance [m]')
    ylabel('Momentum [kg*m/s]')
    F(k) = getframe(h_fig);
end

close(h_fig)
filename = 'gauss_in_quadratic_wigner_repr.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
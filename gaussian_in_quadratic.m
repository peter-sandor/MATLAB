clear;
h_bar = 6.636e-34/2/pi; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
omega = 2*pi*7.5e8;
param_a = mass*omega/(2*h_bar);
alpha_0 = 1e12; 
% If alpha_0 < param_a, wavefunction is the narrowest when it is at the center of
% potential well
% If alpha_0 > param_a, wavefunction is the narrowest when at the extremal
% positions
% If alpha_0 = param_a, then the width of the wavefunction is
% constant in time --> coherent state

% alpha_0 = param_a;
x_0 = 3e-7;
% p_0 = mass*omega*5e-8;
p_0 = 0;
Amp = (2*alpha_0/pi)^(1/4);
time_scale=2*pi/omega; % [s]
N_x = 256;
% N_padded = 2048;
x_rng = 1e-5;
x_stp = x_rng/N_x;
x = ((-x_rng/2+x_stp):x_stp:x_rng/2)'; % [m]
V_quad = 1/2*mass*omega^2*x.^2;
N_delay = 64;
delay_rng = 2*time_scale;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
x_t = x_0*cos(omega*delay)+p_0/(mass*omega)*sin(omega*delay);
p_t = p_0*cos(omega*delay)-mass*omega*x_0*sin(omega*delay);
alpha_t = param_a*(alpha_0*cos(omega*delay)+i*param_a*sin(omega*delay))./(param_a*cos(omega*delay)+i*alpha_0*sin(omega*delay));
gamma_t = (p_t.*x_t-p_0*x_0)/2+i*h_bar/2*log((i*alpha_0*sin(omega*delay)+param_a*cos(omega*delay))/param_a);

scnd_mmnt=zeros([N_delay,1]);
dsprsn=zeros([N_delay,1]);
% Dispersion
h_fig = figure;
for k = 1:length(delay)
% k = 64;
    t=delay(k); % [s]
    Psi = Amp.*exp(-alpha_t(k).*(x-x_t(k)).^2+i/h_bar*p_t(k)*(x-x_t(k))+i/h_bar.*gamma_t(k));
%     Psi_padded = [zeros([floor((N_padded-N_x)/2),1]); Psi; zeros([floor((N_padded-N_x)/2),1])];
%     Phase = mass*alpha_0/sqrt(mass^2+(2*h_bar*alpha_0*t)^2)*sin(2*h_bar*alpha_0*t/mass)*x.^2;
    scnd_mmnt(k) = sqrt(sum(x.^2.*abs(Psi).^2)-sum(x.*abs(Psi).^2).^2);
%     figure;
    subplot(211)
    plot(x,abs(Psi),'k.')
    xlabel('Physical distance [m]')
%     axis([min(x),max(x),-1e4,5e4])
    hold on;
    plot(x,imag(Psi),'b')
    plot(x,real(Psi),'r')
    hold off;
    subplot(212)
    plot(x,V_quad,'g')
    xlabel('Physical distance [m]')
    F(k) = getframe(h_fig);
end
close(h_fig)
filename = 'gaussian_in_quadratic_2';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
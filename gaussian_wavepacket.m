clear;
alpha_0 = 1e18;
h_bar = 6.636e-34/2/pi; % [m^2*kg/s)]
mass = 9.11e-31; % [kg]
Amp = (2*alpha_0/pi)^(1/4);
time_scale=mass/(2*alpha_0*h_bar); % [s]
N_x = 256;
% N_padded = 2048;
x_rng = 1e-7;
x_stp = x_rng/N_x;
x = ((-x_rng/2+x_stp):x_stp:x_rng/2)'; % [m]
N_delay = 64;
delay_rng = 1e-13;
delay_stp = delay_rng/(N_delay);
delay = ((-delay_rng/2+delay_stp):delay_stp:delay_rng/2)'; % [s]
scnd_mmnt=zeros([N_delay,1]);
dsprsn=zeros([N_delay,1]);
% Dispersion
h_fig = figure;
for k = 1:length(delay)
    t=delay(k); % [s]
    alpha_t = alpha_0./(1+2*i*h_bar*alpha_0.*t/mass);
    gamma_t = i*h_bar/2*log(alpha_0./(alpha_t));
    Psi = Amp.*exp(-alpha_t.*x.^2+i/h_bar.*gamma_t);
%     Psi_padded = [zeros([floor((N_padded-N_x)/2),1]); Psi; zeros([floor((N_padded-N_x)/2),1])];
    Phase = mass*alpha_0/sqrt(mass^2+(2*h_bar*alpha_0*t)^2)*sin(2*h_bar*alpha_0*t/mass)*x.^2;
    scnd_mmnt(k) = sqrt(sum(x.^2.*abs(Psi).^2)-sum(x.*abs(Psi).^2).^2);
%     wvnumber_rng = 1/x_stp;
%     wvnumber_stp = 1/x_rng;
    wvnumber_rng = 5/x_stp;
    wvnumber_stp = 5/x_rng;
    wvnumber = ((-wvnumber_rng/2+wvnumber_stp):wvnumber_stp:wvnumber_rng/2)'; % [1/m]
    Psi_spectrum = Amp/(2*sqrt(pi*alpha_0))*exp(-wvnumber.^2/(4*alpha_0)-i*h_bar*wvnumber.^2*t/(2*mass));
    dsprsn(k) = sqrt(sum((h_bar/mass)^2*wvnumber.^2.*abs(Psi_spectrum).^2)-sum(h_bar/mass*wvnumber.*abs(Psi_spectrum).^2).^2); % [m/s];
    prop = 0;
    if k>1
        prop = abs(scnd_mmnt(k)-scnd_mmnt(k-1))/delay_stp/dsprsn(k-1);
    end
    
    subplot(2,2,3)
    plot(wvnumber,abs(Psi_spectrum),'k.')
%     axis([min(wvnumber),max(wvnumber),-3e-6,5e-6])
    xlabel('Wavenumber k, [1/m]')
    hold on
    plot(wvnumber,real(Psi_spectrum),'b')
    plot(wvnumber,imag(Psi_spectrum),'r')
    hold off
    subplot(2,2,4)
    plot(wvnumber,180/pi*h_bar*t*wvnumber.^2/(2*mass),'m')
%     axis([min(wvnumber),max(wvnumber),-2e5,+2e5])
    %plot(wvnumber,180/pi*unwrap(angle(Psi_spectrum))-min(180/pi*unwrap(angle(Psi_spectrum))),'g')
%     hold on;
    
%     figure;
    subplot(2,2,1)
    plot(x,abs(Psi),'k.')
    xlabel('Physical distance [m]')
%     axis([min(x),max(x),-1e4,5e4])
    hold on;
    plot(x,imag(Psi),'b')
    plot(x,real(Psi),'r')
    hold off;
    subplot(2,2,2)
%     %plot(x,180/pi*unwrap(angle(Psi))-min(180/pi*unwrap(angle(Psi))),'g')
%     hold on
    plot(x,180/pi*Phase,'m')
%     axis([min(x),max(x),-12000,+12000])
    title(strcat('Proportionality = ',num2str(prop)));
    % plot(x,180/pi*unwrap(angle(exp(i*Phase))), 'mo')
    F(k) = getframe(h_fig);
end
close(h_fig)
filename = 'gauss_wp.avi';
delete(filename);
movie2avi(F,filename,'compression','Cinepak');
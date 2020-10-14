t_min=0;
t_max=5
t_step=0.01
t=t_min:t_step:t_max;
Topt=5*t_step;
sigma1=0.1;
sigma2=0.1;
delay1=1;
delay2=3;
phase1=0;
phase2=pi;
omega=2*pi/Topt;
field=exp(-(t-delay1).^2/sigma1^2).*cos(omega*t+phase1)+exp(-(t-delay2).^2/sigma2^2).*cos(omega*t+phase2);
intensity=field.^2;

figure;plot(t,intensity)
hold on;plot()
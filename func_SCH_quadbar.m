function FUN = func_SCH_quadbar(delay,P)

mass = 9.1e-31;
h_bar = 6.626e-34;
omega = 2*pi*7.5e8;

FUN = [P(2)/mass; mass*omega^2*P(1);-i*2*h_bar/mass*P(3)^2-i*mass/(2*h_bar)*omega^2;3/2*P(2)^2+h_bar^2/(2*mass)*P(3)+mass*omega^2/2*P(1)^2;];
end
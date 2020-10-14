% simulation for PID control
% Assumptions:
% 1) Linear sequence: We read in Process variable value, calculate the
% Action variable, apply the correction, and then start the cycle over.
% /It may not be exactly the same as how our LabView code works/
% 2) Instantaneous events: We apply correction immediately after we read
% the Process variable value, so there's an immediate reply to external
% drive.
% /It seems a reasonable approximation given that the time between each sampling (delta_t) is ~10 ms, and we see slow drifts on a timescale of 0.1 - 1 s. However, it's an open question what frequencies can the mirror mount mechanics follow./

% The simulation shows that for correcting simple external DC offset, setting the the P parameter only to
% nonzero is enough (e.g.: coeff_prop=0.2-1.0 )
% Also, correcting a slow enough linear drift, the P parameter is sufficient.
% It is also sufficent to correct for slow oscillations; however, only a
% single value of P produces perfect results. (Slow means slow compared to
% sampling rate.)
% With slow oscillations, setting the Integral parameter (coeff_int) to
% nonzero doesn't really help. It just modifies to correct P parameter
% value, but doesn't broaden the range of acceptable P values.
% (E.g.: try coeff_prop=1; coeff_int=0 vs. coeff_prop=0.8; coeff_int=25; )

N = 100;
A = 1;
T_period = 25;
omega = 2*pi/T_period;
procvar_ext = A*sin(omega*(0:100));  % oscillatory external process variable
% procvar_ext = .7*ones([N,1]); % DC external offset case
% procvar_ext = 0:1/(N-1):1; % slow linear change
setpoint = 0.5*ones([N 1]);
coeff_prop = 0.677;
coeff_int = 50;
coeff_der = 0;
T_int = 100;
T_der = 100;
delta_t = 1;

procvar = zeros([N,1]);
procvar(1) = procvar_ext(1);
error = zeros([N,1]);
error(1) = setpoint(1)-procvar(1);
action = zeros([N,1]);
for k=2:N
    procvar(k) = procvar(k-1) + (procvar_ext(k)-procvar_ext(k-1));
    error(k) = setpoint(k)-procvar(k);
    action(k) = coeff_prop*(error(k)+coeff_int/T_int*(error(k)+error(k-1))*delta_t/2+coeff_der*T_der*(error(k-1)-error(k))/delta_t);
    procvar(k) = procvar(k) + action(k);
end

% figure;
plot(setpoint,'k-')
ylim([-1.2 1.2])
hold on;plot(procvar,'r-')
hold on;plot(procvar_ext,'b')
legend('S.P.','P.V. w/ control','P.V. w/o. control')
hold off;
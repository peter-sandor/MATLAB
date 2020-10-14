%% Solve problem of particle in time and space-dependent time-harmonic potential

eps2=1e-5;
z_max = 9e-14;
eps=z_max*1e-3;
z_min = 0;
N = 800;
z_stp = (z_max-z_min)/(N-1);
z = [eps; (z_min+z_stp:z_stp:z_max)'];
param(1)=1e-30; % particle mass
param(2)=2e-9; % force amplitude
param(3)=1e-3; % force scale length
param(4)=2*pi*4e14; % force angular frequency
param(5)=30; % pulse duration

% initial conditions [x(1);x(2);x(3)]
init_cond = [0,0,0];

[zout, yout] = solve_ponderomotive(z,init_cond,param);

% figure('windowstyle','docked');
plot(zout(:,1), yout(:,2:3))
xlim([0,z_max])


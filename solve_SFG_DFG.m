% Solving coupled wave equations for SFG and DFG
% Source: W. R. Boyd Nonlinear Optics (2nd ed. 2003)

z_max = 10;
z_min = 0;
N = 400;
z_stp = (z_max-z_min)/(N-1);
z = (z_min:z_stp:z_max)';
vec_0 = [1;0;0];

[Z, SFG_ode45] = ode45(@func_SFG, z, vec_0);
[Z, DFG_ode45] = ode45(@func_DFG, z, vec_0);

figure('windowstyle','docked')
subplot(211)
plot(Z, abs(SFG_ode45(:,1)), 'k-')
hold on
plot(Z, abs(SFG_ode45(:,2)), 'r-')

subplot(212)
plot(Z, abs(DFG_ode45(:,1)), 'k-')
hold on
plot(Z, abs(DFG_ode45(:,2)), 'r-')
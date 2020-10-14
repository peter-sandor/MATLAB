function [zout, yout] = solve_ponderomotive(zin,init_cond,param)

    function xdot = func_ponderomotive(z,x)
        xdot = [1; x(3)/param(1); param(2)*exp(-x(2)^2/param(3)^2-x(1)^2/param(5)^2)*cos(param(4)*x(1))]; % the defining system of 1st order ODE-s
    end

[zout(:,1), yout(:,:,1)] = ode45(@func_ponderomotive,zin,init_cond); % MATLAB solver
%[zout(:,2), yout(:,:,2)] = RK4(@func_ponderomotive,zin,init_cond); % simple, Home-built 4th order Runge-Kutta solver
end
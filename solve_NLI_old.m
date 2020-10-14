function [Z, NLI_ode45] = solve_NLI(z,init_cond,param)

    function FUN = func_NLI(z,A)

    % Coupled equations for field amplitudes for SFG
    % Source: Boyd W. R. Nonlinear Optics (2nd ed., 2003), p. 73-74.

    K1 = i*param(1);
    K2 = i*param(2);
    K3 = i*param(3);

    FUN = [K1*A(3)*conj(A(2))*exp(-i*param(4)*z); K2*A(3)*conj(A(1))*exp(+i*param(4)*z); K3*A(1)*A(2)*exp(+i*param(4)*z);];
    end

[Z, NLI_ode45] = ode45(@func_NLI,z,init_cond);
end
function FUN = func_SFG(z,A)

% Coupled equations for field amplitudes for SFG
% Source: Boyd W. R. Nonlinear Optics (2nd ed., 2003), p. 79.

K1 = 1i;
K2 = 1i;

FUN = [K1*A(2)*exp(-1i*A(3)*z); K2*A(1)*exp(+1i*A(3)*z); 0];
end
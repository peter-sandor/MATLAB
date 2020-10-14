function FUN = func_DFG(z,A)

% Coupled equations for field amplitudes for DFG
% Source: Boyd W. R. Nonlinear Optics (2nd ed., 2003), p. 84.

K1 = 1i;
K2 = 1i;

FUN = [K1*conj(A(2))*exp(-1i*A(3)*z); K2*conj(A(1))*exp(+1i*A(3)*z); 0];
end
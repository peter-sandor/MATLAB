function FUN = func_NLI(z,A)

% Coupled equations for field amplitudes for SFG
% Source: Boyd W. R. Nonlinear Optics (2nd ed., 2003), p. 73-74.

% c = 2.99792e8;
% omega_1 = ;
% omega_2 = ;
% omega_3 = ;
% k_1 = omega_1/c;
% k_2 = omega_2/c;
% k_3 = omega_3/c;

% neglecting refractive indices, the relative koefficients are set for the
% following scenario:
% 800 nm + 1200 nm --> 480 nm

K1 = i;
K2 = i*1.5;
K3 = i*2.5;

FUN = [K1*A(3)*conj(A(2))*exp(-i*A(4)*z); K2*A(3)*conj(A(1))*exp(+i*A(4)*z); K3*A(1)*A(2)*exp(+i*A(4)*z); 0;];
end
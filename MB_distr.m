function p = MB_distr(E,T)
% Calculate the Maxwell-Boltzmann distribution for the energies 'E' at
% temperature 'T'.
% Use SI units
kB=1.380648e-23; % [J/K]
p=2/sqrt(pi)*sqrt(E).*(1/kB/T).^(3/2).*exp(-E/kB/T);
end
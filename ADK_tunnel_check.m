% Li 2s
Z=3; n_star=2-0.41;
% Na 5s
% Z=3; n_star=5-1.352;
% K 4s
% Z=19; n_star=4-2.23;
% Cs 6s
% Z=55; n_star=6-4.13;

fieldmax=1;
field=0:fieldmax/99:fieldmax;
w=1./field.^(2*n_star-3/2).*exp(-2/3*Z^3/n_star^4./field);
figure;plot(field,w,'k.-');
xlabel('E [atomic units]')
ylabel('Ionization probability')
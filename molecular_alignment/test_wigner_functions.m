% test Wigner functions
J=5;
M=1;
a=wigner_d_P(theta,J);
c=LegendreP(cos(theta),J);
b=wigner_d(theta,[J M 0]);
figure;hold on; 
plot(theta,squeeze(a(:,J+1,J+1+M,J+1)),'k');
plot(theta,b,'r');
plot(theta,c(:,J+1,J+1+M)/max(c(:,J+1,J+1+M)),'b');
legend('Wigner D (my code)','Wigner D (Bob''s code)','Legendre poly')
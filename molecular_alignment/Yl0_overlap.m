maxJ=20;
Jvec=0:maxJ;
Ntheta=200;
theta=map2colvec(0:pi/(Ntheta-1):pi); % array for angle theta
legendre_poly=LegendreP(cos(theta),maxJ);
Yl0=zeros([Ntheta maxJ+1]);
for indJ=1:maxJ+1
    M=0;
    J=Jvec(indJ);
    Yl0(:,indJ)=sqrt((2*J+1)*factorial(J-M)/4/pi/factorial(J+M))*legendre_poly(:,indJ,maxJ+1+M);
end
Nbasis=10;
for indJ=1:maxJ+1
    norm0(indJ)=2*pi*trapz(theta,Yl0(:,indJ).^2.*sin(theta));
    norm2(indJ)=trapz(theta,squeeze(legendre_poly(:,indJ,21)).^2.*sin(theta));
end
str_legend{1}='1';
norm1=2./(2*Jvec+1);%./((2*Jvec+1)/4/pi);
basis = zeros([Ntheta Nbasis]);
for ind1=1:Nbasis
    if ind1>1
        basis(:,ind1)=sin((ind1-1)*theta).^2;
    else
        basis(:,1)=ones([Ntheta 1]);
    end
    str_legend{ind1}=['sin^2(' num2str((ind1-1)) '\theta)'];
    for indJ=1:maxJ+1
        coeff(indJ,ind1)=2*pi*trapz(theta,basis(:,ind1).*Yl0(:,indJ).*sin(theta));
    end
end

figure;plot(Jvec,coeff,'o');
xlabel('J for Y_{J0}(\theta)');
ylabel('expansion coefficients')
legend(str_legend)
%%
indB=5;
calced=zeros([length(Yl0) 1]);
for indJ=1:maxJ+1
    calced = calced + coeff(indJ,indB)*Yl0(:,indJ);
end
figure;plot(theta,[basis(:,indB),calced],'-')
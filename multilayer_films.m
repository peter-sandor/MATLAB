% calculate reflection and transmission coefficients for thin multilayer
% films
% microns and picoseconds are used for units of distance
c=300; % [um/ps]
lambda=0.780;
theta0=20; % in degrees
theta0_rad=theta0/180*pi;
n0=1;
nS=1.5;
thetaS_rad=asin(n0*theta0_rad/nS);
polarization=1; % 0 for p, 1 for s
layers=[2 1.5; 3 1.1]; % [thickness_1 refractive_index_1; thickness_2 refractive_index_2; ...]

Mtot=eye(2);
for ind1=1:size(layers,1)
    theta_rad(ind1)=asin(n0*theta0_rad/layers(ind1,2));
    delta(ind1)=2*pi*layers(ind1,1)*layers(ind1,2)*cos(theta_rad(ind1));
    if polarization==0
        gamma(ind1)=layers(ind1,2)/c/cos(theta_rad(ind1));
    else 
        gamma(ind1)=layers(ind1,2)/c*cos(theta_rad(ind1));
    end
    M(:,:,ind1)=[cos(delta(ind1)),i*sin(delta(ind1))/gamma(ind1);i*gamma(ind1)*sin(delta(ind1)),cos(delta(ind1))];
    Mtot=Mtot*M(:,:,ind1);
end

if polarization==0
    term1=cos(thetaS_rad*Mtot(1,1)+nS/c*Mtot(1,2));
    term2=cos(thetaS_rad*Mtot(2,1)+nS/c*Mtot(2,2));
    refl=(term2*cos(theta0_rad)-n0/c*term1)/(term2*cos(theta0_rad)+n0/c*term1);
    transm=2*n0*cos(theta0_rad)/(term1*cos(theta0)+n0/c*term1);
else
    gamma0=n0/c*cos(theta0_rad);
    gammaS=nS/c*cos(thetaS_rad);
    refl=(gamma0*Mtot(1,1)+gamma0*gammaS*Mtot(1,2)-Mtot(2,1)-Mtot(2,2)*gammaS)/(gamma0*Mtot(1,1)+gamma0*gammaS*Mtot(1,2)+Mtot(2,1)+Mtot(2,2)*gammaS);
    transm=2*gamma0/(gamma0*Mtot(1,1)+gamma0*gammaS*Mtot(1,2)+Mtot(2,1)+Mtot(2,2)*gammaS);
end


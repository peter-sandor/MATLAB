N=50;
l=1:10;
coeffs=ones(size(l));
% coeffs=[1 0 1 0 1 0 1 0 1 0 1];
m=0;
theta=0:pi/(N-1):pi;
phi=0:2*pi/(N-1):2*pi;
r=ones([1 N]);
for ind3=1:length(l)
    Ylm(:,:,ind3)=sph_harmonics(l(ind3),m,theta,phi);
end
for ind1=1:N
    for ind2=1:N
        [x(ind1,ind2),y(ind1,ind2),z(ind1,ind2)] = sph2cart(phi(ind2),theta(ind1),1);
    end
end
figure;
% imagesc(phi/pi*180,theta/pi*180,double(real(Ylm)>=0))
imagesc(phi/pi*180,theta/pi*180,abs(sum(Ylm.*permute(extend(extend(coeffs,N),N),[2 3 1]),3)))
% surf(x,y,z,double(real(Ylm)>=0));
% colormap gray
% shading flat
%%
N=50;
l=0:10;
m=l;%0*ones(size(l));
% coeffs=ones(size(l));
coeffs=[1 0 1 0 1 0 1 0 1 0 1];
[x2,y2,z2] = sphere(N);
for ind1=1:N
    for ind2=1:N
        [phi2(ind1,ind2),theta2(ind1,ind2),r]=cart2sph(x2(ind1,ind2),y2(ind1,ind2),z2(ind1,ind2));
        theta2(ind1,ind2)=theta2(ind1,ind2)+pi/2;
        for ind3=1:length(l)
            Ylm2(ind1,ind2,ind3)=sph_harmonics(l(ind3),m(ind3),theta2(ind1,ind2),phi2(ind1,ind2));
        end
    end
end
% figure;
% surf(x2,y2,z2,double(real(Ylm2)>0))
surf(x2,y2,z2,abs(sum(Ylm2.*permute(extend(extend(coeffs,N),N),[2 3 1]),3)).^2)
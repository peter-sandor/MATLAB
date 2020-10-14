% rotational constants
B=0.203; % [1/cm]
D=3.46e-8; % [1/cm]
T=10;

maxJ=20;
J=0:maxJ;
h=6.626e-34; % [J*s]
c=3e8; % [m/s]
E=h*c*100*(B*J.*(J+1)+D*J.^2.*(J+1).^2);
p=MB_distr(E,T);
coeffsJ=p/sum(p);

figure;
subplot(121)
plot(J,E/1.6e-19)
subplot(122)
bar(J,coeffsJ)

%%
% N=50;
% m=0;
% theta=0:pi/(N-1):pi;
% phi=0:2*pi/(N-1):2*pi;
% r=ones([1 N]);
% for ind3=1:length(J)
%     Ylm(:,:,ind3)=sph_harmonics(J(ind3),m,theta,phi);
% end
% for ind1=1:N
%     for ind2=1:N
%         [x(ind1,ind2),y(ind1,ind2),z(ind1,ind2)] = sph2cart(phi(ind2),theta(ind1),1);
%     end
% end
% figure;
% imagesc(phi/pi*180,theta/pi*180,double(real(Ylm)>=0))
% imagesc(phi/pi*180,theta/pi*180,abs(sum(Ylm.*permute(extend(extend(coeffsJ,N),N),[2 3 1]),3)))
% surf(x,y,z,double(real(Ylm)>=0));
% colormap gray
% shading flat
%%
N=50;
[x2,y2,z2] = sphere(N);
Ylm=zeros([N N maxJ+1 2*maxJ+1]);
for ind1=1:N
    for ind2=1:N
        [phi2(ind1,ind2),theta2(ind1,ind2),r]=cart2sph(x2(ind1,ind2),y2(ind1,ind2),z2(ind1,ind2));
        theta2(ind1,ind2)=theta2(ind1,ind2)+pi/2;
        for ind3=1:length(J)
            m=-J(ind3):J(ind3);
            for ind4=1:length(m)
                Ylm(ind1,ind2,ind3,ind4)=sph_harmonics(J(ind3),m(ind4),theta2(ind1,ind2),phi2(ind1,ind2));
            end
        end
    end
    disp([num2str(ind1) '/' num2str(N)]);
end
%%
coeffs=zeros([maxJ+1 2*maxJ+1]);
% for ind3=1:maxJ+1
%     coeffs(ind5,ind5)=coeffsJ(ind5);
% end
delay=(0:1:90)*1e-12;
Ndelay=length(delay);
for ind5=1:Ndelay
    for ind3=1:length(J)
        for ind4=1:(2*J(ind3)+1)
            coeffs(ind3,ind4,ind5)=mod(J(ind3),2)*coeffsJ(ind3)/(2*J(ind3)+1)*exp(1i*1*rand)*exp(-1i*2*pi*E(ind3)/h*delay(ind5)); % weighing different 'm'-s equally, but with random phase
        end
    end
end
Exp_val_cos2=zeros([Ndelay 1]);
Ylm_sum=zeros([N N Ndelay]);
for ind5=1:Ndelay
    Ylm_sum(:,:,ind5)=sum(sum(Ylm.*permute(extend(extend(coeffs(:,:,ind5),N),N),[3 4 1 2]),3),4);
    Exp_val_cos2(ind5)=abs(sum(sum(Ylm_sum(:,:,ind5).*sin(theta2).*cos(theta2).^2.*conj(Ylm_sum(:,:,ind5))))).^2/abs(sum(sum(Ylm_sum(:,:,ind5).*sin(theta2).*conj(Ylm_sum(:,:,ind5))))).^2;
    Exp_val_cos4(ind5)=abs(sum(sum(Ylm_sum(:,:,ind5).*sin(theta2).*cos(theta2).^4.*conj(Ylm_sum(:,:,ind5))))).^2/abs(sum(sum(Ylm_sum(:,:,ind5).*sin(theta2).*conj(Ylm_sum(:,:,ind5))))).^2;
    disp([num2str(ind5) '/' num2str(Ndelay)]);
end

figure;
% surf(x2,y2,z2,double(real(Ylm2)>0))
surf(x2,y2,z2,abs(Ylm_sum(:,:,1)).^2);
shading flat;

figure;plot(delay,Exp_val_cos2)
%% evaluate <cos(theta)^2> and <cos(theta)^4>
% 
% Exp_val_cos2=abs(sum(sum(Ylm_sum.*sin(theta2).*cos(theta2).^2.*conj(Ylm_sum)))).^2/abs(sum(sum(Ylm_sum.*sin(theta2).*conj(Ylm_sum)))).^2;
% Exp_val_cos4=abs(sum(sum(Ylm_sum.*sin(theta2).*cos(theta2).^4.*conj(Ylm_sum)))).^2/abs(sum(sum(Ylm_sum.*sin(theta2).*conj(Ylm_sum)))).^2;
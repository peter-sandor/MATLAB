function draw_sph_grid(imgin,center,Ntheta,Nrad)

N=size(imgin);
% center=[151 151];
cutoff=0.003;
% Nbin=20;
% Nsector=20;
dist=sort([center(1) N(2)-center(1) center(2) N(1)-center(2)],'descend');
Rmax=dist(ceil(N(1)*cutoff));
flag1=max(cumprod(size(Nrad)));
if flag1==1
    Rbins=map2colvec(0:Rmax/(Nrad-1):Rmax);
else
    Rbins=Nrad;
end
Sbins=(0:2*pi/Ntheta:2*pi*(1-1/(Ntheta)))+pi/2;
imagesc(imgin)
for ind1=1:Ntheta
    if abs(cos(Sbins(ind1)))>=1e-4
        xvec=center(2):cos(Sbins(ind1))*Rmax:cos(Sbins(ind1))*Rmax+center(2);
    else
        xvec=[center(2) center(2)];
    end
    if abs(sin(Sbins(ind1)))>=1e-4
        yvec=center(1):sin(Sbins(ind1))*Rmax:sin(Sbins(ind1))*Rmax+center(1);
    else
        yvec=[center(1) center(1)];
    end
    line(xvec,yvec,'color','w')
%     disp(num2str(ind1))
end
M=50;
if flag1==1
    for ind1=1:length(Rbins)
        xvec=-Rbins(ind1)+center(2):2*Rbins(ind1)/(M-1):Rbins(ind1)+center(2);
        yvec=sqrt((Rbins(ind1))^2-(xvec-center(2)).^2)+center(1);
        line(xvec,yvec,'color','w','linewidth',2)
        line(xvec,2*center(1)-yvec,'color','w','linewidth',2)
    %     disp(num2str(ind1))
    end
else
    load COLORS_line;
    for ind1=1:size(Rbins,1)
        for ind2=1:size(Rbins,2)
            xvec=-Rbins(ind1,ind2)+center(2):2*Rbins(ind1,ind2)/(M-1):Rbins(ind1,ind2)+center(2);
            yvec=sqrt((Rbins(ind1,ind2))^2-(xvec-center(2)).^2)+center(1);
            line(xvec,yvec,'color',COLORS_line(ind1,:),'linewidth',2)
            line(xvec,2*center(1)-yvec,'color',COLORS_line(ind1,:),'linewidth',2)
        %     disp(num2str(ind1))
        end
    end
end
end
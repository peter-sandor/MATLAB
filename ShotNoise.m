function imgout = ShotNoise(imgin,sigma,level,N)
M=size(imgin);
imgin=imgin/max(max(imgin));
[X,Y]=meshgrid(1:M(2),1:M(1));
amp=level*max(max(imgin));
imgout=zeros(size(imgin));
ind1=0;
while ind1<=N
    newcoord(1)=M(1)*rand;
    if newcoord(1)<1
        newcoord(1)=1;
    end
    newcoord(2)=M(2)*rand;
    if newcoord(2)<1
        newcoord(2)=1;
    end    
    flrd=floor(newcoord);
    rmndr=mod(newcoord,1);
    value=interpP(imgin(flrd(1):flrd(1)+1,flrd(2):flrd(2)+1),rmndr);
    if value>=rand
        shot=amp*exp(-(X-newcoord(2)).^2/sigma-(Y-newcoord(1)).^2/sigma);
        imgout=imgout+shot;
        ind1=ind1+1;
    end
end
end
function synth_pic = synth_image(coords,picsize,sigma)

Ny=picsize(1);
Nx=picsize(2);
vecx=1:Nx;
vecy=1:Ny;
[X,Y]=meshgrid(1:Nx,1:Ny);
synth_pic=zeros(picsize);
for ind1=1:size(coords,1)
    synth_pic = synth_pic + 1/sigma^2/pi*exp(-((X-coords(ind1,1)).^2+(Y-coords(ind1,2)).^2)/sigma^2);
end
end
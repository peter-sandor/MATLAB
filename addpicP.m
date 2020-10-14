function out = addpicP(extension,main,bckg)
% extension is file extension
% main is an array with the indices of main images
% bckg is an array with the indices of background images

files = dir(strcat('*.',extension));
for ind = main
    pic = monopic(imread(files(ind).name));
    if ind == main(1)
        sumpic = pic;
    else
        sumpic = sumpic + pic;
    end;
end;
sumpic = sumpic /length(main);
save sumpic.mat sumpic;

for ind = bckg
    pic = monopic(imread(files(ind).name));
    if ind == bckg(1)
        bkgpic = pic;
    else
        bkgpic = sumpic + pic;
    end;
end;

bkgpic = bkgpic /(length(bckg));
save bkgpic.mat bkgpic
load detD
out = (sumpic - bkgpic)./detD;



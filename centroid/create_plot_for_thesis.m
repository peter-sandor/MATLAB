% for thesis

pic=monopic(imread('M:\data2015\2015_02_05_R\test\pic_14.jpg'));
pic2=pic(400:900,900:1600);
params.thrs=33;
params.radius=2;
params.minlistlength=4;
params.filtersize=60;
params.max_coord_diff=30;
thrsd=pic2;
thrsd(thrsd<37)=0;
coords2=centroid_M(pic2,params);
hitimg=create_hitimg(circshift(coords2,[0 1]),4,size(pic2));

figure;
subplot(121)
imagesc(pic2)
title('raw image')
% subplot(312)
% imagesc(thrsd)
subplot(122)
imagesc(hitimg)
title('synthesized image')
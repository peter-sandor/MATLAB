pic=monopic(imread('M:\2015_02_05_R\test\pic_14.jpg'));
% pic=monopic(imread('M:\data2015\2015_03_30_R\for_code_testing3\pic1.png'));
% pic=monopic(imread('M:\data2015\2015_03_30_R\for_code_testing2\pic_6.jpg'));
% pic=monopic(imread('testimages\sub_pic_15.jpg'));

% tic;
% coords=centroid_subimage(pic,32,[1 1]);
% t=toc;

params.thrs=33;
params.radius=2;
params.minlistlength=4;
params.filtersize=60;
params.max_coord_diff=30;

% tic;
% coords=centroid_basic(pic,params);
% t1=toc;

tic;
coords2=centroid_M(pic,params);
t2=toc;

% hitimg=create_hitimg(coords,4,size(pic));
% dfigure;
% subplot(121)
imagesc(pic)
hold on
plot(coords(:,1),coords(:,2),'r+')
plot(coords2(:,1),coords2(:,2),'ko')
hold off
% subplot(122)
% imagesc(hitimg)

%%
pic=double(color2gray(imread('testimages/sub_pic_9.jpg')));
% pic=color2gray(imread('testimages/pic113.jpg'));
% pic=color2gray(imread('testimages/just_a_hit2.png'));
maxval=double(max(max(pic)));
hist1=hist(reshape(double(pic),[size(pic,1)*size(pic,2) 1]),maxval);
range=1:maxval;
meanval=sum(range.*hist1)/sum(hist1);
FIT1=ezfit(range,hist1,['y=B*A^x/factorial(x)*exp(-A);A=' num2str(meanval) ';B=' num2str(max(hist1)) ';']);
fitted=FIT1.m(2)*FIT1.m(1).^range./factorial(range)*exp(-FIT1.m(1));
figure;bar(range,hist1);
hold on;plot(range,fitted,'r');
tic;
coords = centroid_bwconn(pic,2.5*FIT1.m(1));
t3=toc;
pic_synth=synth_image(coords,size(pic),4);
figure;
subplot(121)
imagesc(pic)
hold on;plot(coords(:,1),coords(:,2),'r*')
subplot(122)
imagesc(pic_synth)
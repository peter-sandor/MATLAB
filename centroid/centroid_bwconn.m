function coords = centroid_bwconn(pic,thrs)
% pic=color2gray(imread('testimages/3hits.jpg'));
% thrs=50;
% pic=pic;
pic(pic<thrs)=0;
pic(pic>=thrs)=1;
CC = bwconncomp(pic);%Finds objects
stats=regionprops(CC,'centroid');
N=length(stats);
for ind1=1:N
    coords(ind1,:)=stats(ind1).Centroid;
end
% figure;imagesc(pic)
% hold on;plot(coords(:,1),coords(:,2),'r*')
end
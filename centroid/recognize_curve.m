function [x,y]=recognize_curve(pic,thrs)

% params.thrs=round((max(max(pic))-min(min(pic)))/2);
params.thrs=thrs;
params.radius=10;
params.minlistlength=10;
params.filtersize=60;
params.max_coord_diff=25;
params.swap_xy=1;

% pic=monopic(imread('curve_for_test.jpg'));
pic=flipdim(pic,1);
[a,b]=centroid_M(pic,params);
S=[];
for ind1=1:length(b)
    S=cat(1,S,b{ind1});
end
[N1 N2]=size(pic);
picr=zeros([N1 N2]);
for ind1=1:size(S,1)
    picr(S(ind1,2),S(ind1,1))=S(ind1,3);
end

S2=S;
S2(:,3)=[];
S2=sort(S2,1,'ascend');
[x,y]=consolidator(S(:,1),S(:,2),'mean');
end
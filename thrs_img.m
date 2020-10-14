function varargout = thrs_img(img,thrs,radius,minlistlength,filtersize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N=size(img);
ind3=1;
pts=zeros([10000 3]);
for ind1=2:N(1)-1
    for ind2=2:N(2)-1
%         if (img(ind1,ind2)+img(ind1-1,ind2)+img(ind1,ind2-1)+img(ind1+1,ind2)+img(ind1,ind2+1))>5*thrs
         if img(ind1,ind2)>thrs
            pts(ind3,:)=[ind1 ind2 img(ind1,ind2)];
            ind3=ind3+1;
        end
    end
end
pts(pts(:,1)==0,:)=[];
if size(pts,1)>1
    pts2=[];
    ind2=1;
    Npts=size(pts,1);
    for ind1=1:Npts
        range_min=max(1,ind1-filtersize);
        range_max=min(Npts,ind1+filtersize);
        temp=double(((pts(ind1,1)-pts(range_min:range_max,1)).^2+(pts(ind1,2)-pts(range_min:range_max,2)).^2)<=radius^2);
        if sum(temp)>minlistlength
            pts2(ind2,:)=pts(ind1,:);
            ind2=ind2+1;
        end
    end
else
    pts2=[];
end
varargout{1}=pts2;
end
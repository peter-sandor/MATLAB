function imageout = foldimage(imagein,center)
% this code takes an image and the coordinates of its center, then folds it
% into one quadrant. Useful for processing photoelectron or photoion VMI
% images.
% center: [x,y] a.k.a. [column,row] (can be non-integer)

M=flipdim(size(imagein),2);
if (round(center(1))~=center(1)) || (round(center(2))~=center(2)) % shift and resample image if center coordinates are not integers
    [X,Y]=meshgrid(1:M(1),1:M(2));
    [XI,YI]=meshgrid((1+mod(center(1),1)):1:(M(1)-1+mod(center(1),1)),(1+mod(center(2),1)):1:(M(2)-1+mod(center(2),1)));
    imagein = interp2(X,Y,imagein,XI,YI,'linear');
    center=[center(1)-mod(center(1),1) center(2)-mod(center(2),1)];
end
M=flipdim(size(imagein),2); % update image size
a=[center(1) M(1)-center(1)+1 center(2) M(2)-center(2)+1]; %[xl xr yt yb]
[b,c]=sort(a(1:2),2,'descend');
if c(1)==1
    permvec(1)=1;
else
    permvec(1)=2;
end
[b,c]=sort(a(3:4));
if c(1)==1
    permvec(2)=3;
else
    permvec(2)=4;
end
permvec=sort(permvec(1:2));
if findvec(permvec,[1 3]) % [xl yt] --> first quadrant
    center2=center+[-1 -1];
    N=min([center2(1) center2(2)]); % N = size of one quadrant
elseif findvec(permvec,[2 3])  % [xr yt] --> second quadrant
    center2=center+[+1 -1];
    N=min([M(1)-center2(1)+1 center2(2)]);
elseif findvec(permvec,[1 4])  % [xl yb] --> third quadrant
    center2=center+[-1 +1];
    N=min([center2(1) M(2)-center2(2)+1]);
elseif findvec(permvec,[2 4])  % [xr yb] --> fourth quadrant
    center2=center+[+1 +1];
%     N=10*floor(min([M(1)-center2(1)+1 M(2)-center2(2)+1])/10);
    N=min([M(1)-center2(1)+1 M(2)-center2(2)+1]);
%     N2=10*(max([M(1)-center2(1)+1 M(2)-center2(2)+1])/10);
end

cropped=imagein(center(2)-N:center(2)+N,center(1)-N:center(1)+N); % crop the image
fold1=zeros([N 2*N+1]);
% fold2=zeros([N2 2*N2+1]);
fold1=cropped(1:N+1,:)+flipdim(cropped(N+1:2*N+1,:),1);
% fold2=imagein(1:center(2),:)+flipdim(cropped(center(2)+1:2*N+1,:),1);
imageout=zeros(N);
imageout=flipdim(flipdim(fold1(:,1:N+1)+flipdim(squeeze(fold1(:,N+1:2*N+1)),2),1),2);
end
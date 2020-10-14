function imageout = foldimage(imagein,center,quads)
% this code takes an image and the coordinates of its center, then folds it
% into one quadrant. Useful for processing photoelectron or photoion VMI
% images.
% center: [x,y] a.k.a. [column,row] (can be non-integer)
% quadrant convention: top left: 1, top right: 2, bottom left: 3, bottom
% right: 4

M=flipdim(size(imagein),2);
if (round(center(1))~=center(1)) || (round(center(2))~=center(2)) % shift and resample image if center coordinates are not integers
    [X,Y]=meshgrid(1:M(1),1:M(2));
    [XI,YI]=meshgrid((1+mod(center(1),1)):1:(M(1)-1+mod(center(1),1)),(1+mod(center(2),1)):1:(M(2)-1+mod(center(2),1)));
    imagein = interp2(X,Y,imagein,XI,YI,'linear');
    center=[center(1)-mod(center(1),1) center(2)-mod(center(2),1)];
end
M=flipdim(size(imagein),2); % update image size

quadsizes=[center(1) center(2);
            M(1)-center(1)+1 center(2);
            center(1) M(2)-center(2)+1;
            M(1)-center(1)+1 M(2)-center(2)+1];
quadsizes=quadsizes(quads,:);
N=min(min(quadsizes));
imageout=zeros(N);
for ind1=1:length(quads)
    switch quads(ind1)
        case 1
            imageout=imageout+imagein(center(2)-(N-1):center(2),center(1)-(N-1):center(1));
        case 2
            imageout=imageout+flipdim(imagein(center(2)-(N-1):center(2),center(1):center(1)+(N-1)),2);
        case 3
            imageout=imageout+flipdim(imagein(center(2):center(2)+(N-1),center(1)-(N-1):center(1)),1);
        case 4
            imageout=imageout+flipdim(flipdim(imagein(center(2):center(2)+(N-1),center(1):center(1)+(N-1)),2),1);
    end
end
imageout=flipdim(flipdim(imageout,2),1); % put center to [1 1]
end
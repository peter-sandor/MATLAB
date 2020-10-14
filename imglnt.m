function [lineout,plf] = imglnt(img,A)

% This function calculates the lineout of an image between two specified
% points: A=(x1,y1) and B=(x2,y2). The basic idea is to crop the image to only
% include the rectangle between the above two points, and divide it into N
% strips perpendicular to the line connecting A to B. Here N is the number
% of points in the lineout. Then for each point of the lineout, the value
% is calculated as the weighted average of the values at the points consisting the
% corresponding strip. The weight for each point is chosen to depend
% on the distance between the point and the line connecting A and B; for
% points further apart from the line, the weight is smaller.

A=round(A); % A contains the indices of the two points: A=[x1 y1;x2 y2];
temp=A;
if A(1,1)==A(2,1)
    lineout=map2colvec(img(min(A(1,2),A(2,2)):max(A(1,2),A(2,2)),A(1,1)));
    plf=1;
    if A(2,2)<A(1,2)
        lineout=flipdim(lineout,1);
    end
elseif A(1,2)==A(2,2)
    lineout=map2colvec(img(A(1,2),min(A(1,1),A(2,1)):max(A(1,1),A(2,1))));
    plf=1;
    if A(2,1)<A(1,1)
        lineout=flipdim(lineout,1);
    end
else
    index1_start=uint16(round(min(temp(:,2))));
    index1_end=uint16(round(max(temp(:,2))));
    index2_start=uint16(round(min(temp(:,1))));
    index2_end=uint16(round(max(temp(:,1))));
    if index1_start < 1
        index1_start=1;
    end
    if index2_start < 1
        index2_start=1;
    end
    if index1_end > size(img,1)
        index1_end=size(img,1);
    end
    if index2_end > size(img,2)
        index2_end=size(img,2);
    end 
    img=img(index1_start:index1_end,index2_start:index2_end);
    [M1 M2]=size(img);
    N=max(abs(A(1,1)-A(2,1))+1,abs(A(1,2)-A(2,2))+1);
    lineout=zeros([N 1]);
    plf=sqrt((A(1,1)-A(2,1))^2+(A(1,2)-A(2,2))^2)/max(abs(A(1,1)-A(2,1)),abs(A(1,2)-A(2,2))); % physical length factor
% M1=abs(indices(1,1)-indices(2,1));
% M2=abs(indices(1,2)-indices(2,2));
    m(1)=-(A(2,2)-A(1,2)+1)/(A(2,1)-A(1,1)+1);
    b=0;
% b=(A(2,2)+A(1,2))/(A(2,1)+A(1,1))/2/m(1);
    m(2)=-1/m(1);
    [X,Y] = meshgrid((0:size(img,2)-1),(0:size(img,1)-1));
    Y=-Y;
    % Coordinates are generated from indices based on the displacement
    % vector input A defines: region I: x>0 & y>0; region II: x<0 & y>0;
    % region III: x<0 & y<0; region IV: x>0 & y<0;
	if diff(A(:,2))<0
        Y=Y-min(min(Y)); % Regions III and IV
    end
    if diff(A(:,1))<0 % Regions II and III
        X=X-max(max(X));
    end
    B=Y-m(2)*X;
    sectn_x=(B-b)/(m(1)-m(2));
    sectn_y=1/2*((m(1)+m(2))*sectn_x+b+B);
    dist=sqrt((X-sectn_x).^2+(Y-sectn_y).^2);
    sectn_x=abs(sectn_x);
    sctnx_values=unique(sectn_x);
    sectn_step=(max(sctnx_values)-min(sctnx_values))/(N-1);
    sectn_binned=map2colvec(min(sctnx_values):sectn_step:max(sctnx_values));
    alpha=1;
    lineout=zeros([N 1]);
    for ind1=1:N
        indices1=(sectn_x>=(sectn_binned(ind1)-sectn_step/2));
        indices2=(sectn_x<=(sectn_binned(ind1)+sectn_step/2));
        indices3=logical(indices1.*indices2);
        weights=exp(-alpha*dist(indices3));
        norm=sum(weights);
        lineout(ind1)=sum(img(indices3).*weights)/norm;
    end
end
end
function centerout = getmean(varargin)
% Determine mean of a distribution (pic) within a specific domain (defined
% by A). 'A' has the format: [x1 y1; x2 y2], where the points correspond to
% the two opposite corners of the rectangular domain.
% centerout = [x y] = [column row]

if nargin==0
    pic=get(get(gca,'Children'),'cdata');
    A=ginput(2);
elseif nargin==1;
    pic=varargin{1};
    A=ginput(2);
elseif nargin==2;
    pic=varargin{1};
    A=varargin{2};
end
A=sortrows(A);
A(1,1)=floor(A(1,1));
A(2,1)=ceil(A(2,1));
if A(2,2)>A(1,2)
    A(2,2)=ceil(A(2,2));
    A(1,2)=floor(A(1,2));
    x=A(1,1):A(2,1);
    y=A(1,2):A(2,2);
else A(2,2)=floor(A(2,2));
    A(1,2)=ceil(A(1,2));
    x=A(1,1):A(2,1);
    y=A(2,2):A(1,2);
end
y=map2colvec(y);
pic2=pic(y,x);
centerout=[sum(sum(pic2,1).^2.*x)./sum(sum(pic2,1).^2) sum(sum(pic2,2).^2.*y)./sum(sum(pic2,2).^2)];
centerout=flipdim(centerout,2);
end
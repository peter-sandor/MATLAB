function varargout = plot_FragImgP_Ver4(varargin)

if nargin==1
    ImgSum = zeros(360);
    ImgList=varargin{1};
elseif nargin==2
    ImgList=varargin{1};
    ImgSum = varargin{2};
end

ImgNum = size(ImgList,1);
for j = 1:ImgNum

%         Img = loadbinaryimage2([FileList{k,3} '\pic' num2str(ImgList(j,1)) '.bin']);
%         Img(1,4)=0;
%         Img(1,8)=0;
%         Img = Img.*(Img > 0.01);
    Img = zeros(360);
    cx = ImgList(j,2);
    cy = ImgList(j,1);
    for y = max((round(cy)-2),1):min((round(cy)+2),360)
        for x = max((round(cx)-2),1):min((round(cx)+2),360)
            if x==0 || y==0 || x==360 || y==360
%                     disp([num2str(x) num2str(y)])
                break
            end
            Img(y,x) = exp(-(x-cx)^2/2-(y-cy)^2/2);
        end
    end
    if size(ImgSum,1)~=size(Img,1)
        disp(['j=' num2str(j)]);
    end
    ImgSum = ImgSum + Img;
end
varargout{1}=ImgSum;
end

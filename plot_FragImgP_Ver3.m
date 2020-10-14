function varargout = plot_FragImgP_Ver3(varargin)

if nargin==0
    ImgSum = zeros(360);
    reply='y';
    k=1;
    while reply=='y';
        [filename,path] = uigetfile('*.txt');
 %     FileList{k,3} = uigetdir;
        reply = input([DoubleBackSlash(path) filename ' added. More? Y/N [Y]: '], 's');
        FileList{k}=[path filename];
        k=k+1;
    end
elseif nargin==1
    FileList=varargin{1};
    ImgSum = zeros(360);
elseif nargin==2
    FileList=varargin{1};
    ImgSum = varargin{2};
end

FileNum = length(FileList);
for k=1:FileNum
    ImgList = dlmread(FileList{k});
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
end
varargout{1}=ImgSum;
end

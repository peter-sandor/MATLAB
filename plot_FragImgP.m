function varargout = plot_FragImgP(coords,sigma,imgsize)

% sigma=4;
ImgSum = zeros(imgsize);
N = size(coords,1);
for k=1:N
    Img = zeros(imgsize);
    cx = coords(k,1);
    cy = coords(k,2);
    for y = max((round(cy)-2*sigma),1):min((round(cy)+2*sigma),imgsize(1))
        for x = max((round(cx)-2*sigma),1):min((round(cx)+2*sigma),imgsize(2))
            if x~=0 && y~=0 && x~=imgsize(2) && y~=imgsize(1)
                Img(y,x) = exp(-(x-cx)^2/sigma^2-(y-cy)^2/sigma^2);
            end
        end
    end
    ImgSum = ImgSum + Img;
end
varargout{1}=ImgSum;
end
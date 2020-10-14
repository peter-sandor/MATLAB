function img_out = imagegauss(imsize,center,radius,sigma)

% input parameters:
%   imsize: [sizey sizex] = [rowlength columnlength]
%   center: [y x] = [row column]
%   radius: [inner outer] in pixels
%   sigma: in pixels

img_out=zeros(imsize);
for ind1=1:imsize(1)
	for ind2=1:imsize(2)
        rad=sqrt((ind1-center(1))^2+(ind2-center(2))^2);
        if (rad<=radius(2)) && (rad>=radius(1))
            img_out(ind2,ind1)=1;
        elseif rad>radius(2)
            img_out(ind2,ind1)=exp(-(rad-radius(2))^2/sigma^2);
        elseif rad<radius(1)
            img_out(ind2,ind1)=exp(-(rad-radius(1))^2/sigma^2);
        end
    end
end
img_out=double(img_out);
end
function index_out = imagecirc(imsize,center,radius)

% input parameters:
%   imsize: [sizey sizex] = [columnlength rowlength]
%   center: [y x] = [row column]
%   radius: [inner outer] in pixels

index_out=zeros(imsize);
for ind1=1:imsize(1)
	for ind2=1:imsize(2)
        if (((ind1-center(1))^2+(ind2-center(2))^2)<=radius(2)^2) && (((ind1-center(1))^2+(ind2-center(2))^2)>=radius(1)^2)
            index_out(ind1,ind2)=1;
        end
    end
end
index_out=logical(index_out);
end
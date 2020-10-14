function pic = color2gray(pic0)
% Converts an NxMx3 color image (e.g. from a .bmp) to a grayscale image (NxM double precision)
if size(pic0,3)==3
    pic = squeeze(0.3*double(pic0(:,:,1)) + 0.59*double(pic0(:,:,2)) + 0.11*double(pic0(:,:,3)));
else
    pic = pic0;
end
end
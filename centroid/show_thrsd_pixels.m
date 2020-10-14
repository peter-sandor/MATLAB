function img_out = show_thrsd_pixels(pts,framesize)

img_out=zeros(framesize);
for ind1=1:size(pts,1)
    img_out(pts(ind1,2),pts(ind1,1))=1;%pts(ind1,3);
end
end
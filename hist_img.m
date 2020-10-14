function varargout = hist_img(pic)

pic_array=reshape(pic,[size(pic,1)*size(pic,2) 1]);
pic_hist=hist(pic_array,1:max(pic_array));
figure;bar(pic_hist)
varargout{1}=pic_hist;
end
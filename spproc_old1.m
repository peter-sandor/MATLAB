function spproc
% data: Nx2 array; first column: wavelengths, second column: corresponding
% spectral densities
S=uiimport('-file');
data=S.data;
clear S
hndl1=figure;plot(data(:,2),'k')
title('Choose start and end of background')
temp=ginput(2);
bckg_start=uint16(round(min(temp(1,1),temp(2,1))));
bckg_end=uint16(round(max(temp(1,1),temp(2,1))));
% bckg_start=min(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
% bckg_end=max(uint16(round((temp(:,1)-data(1,1))/(data(size(data,1),1)-data(1,1))*size(data,1))));
if bckg_start < 1
    bckg_start=1;
end
if bckg_end > size(data,1)
    bckg_end=size(data,1);
end        
title('Select region of interest')
temp=ginput(2);
crop_start=uint16(round(min(temp(1,1),temp(2,1))));
crop_end=uint16(round(max(temp(1,1),temp(2,1))));
spectrum=(data(crop_start:crop_end,2)-mean(data(bckg_start:bckg_end,2)))/max(data(crop_start:crop_end,2)-mean(data(bckg_start:bckg_end,2)));
wl=data(crop_start:crop_end,1);
meanwl=sum(wl.*spectrum)/sum(spectrum)
std_dev=std(spectrum)*(max(wl)-min(wl))/length(wl)
plot(wl,spectrum,'k.-')
hold on;plot(wl,0.5*ones(size(spectrum)),'r')
hold off;
title('Choose points for FWHM')
temp=ginput(2);
close(hndl1);
FWHM=abs(diff(temp(:,1)))

end
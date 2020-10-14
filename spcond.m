function varargout = spcond(spectrum)

hndl1=figure;plot(spectrum,'k')
% Set indices for background subtraction
title('Choose start and end of background')
temp=ginput(2);
bkg(1)=uint16(round(min(temp(1,1),temp(2,1))));
bkg(2)=uint16(round(max(temp(1,1),temp(2,1))));
if bkg(1) < 1
    bckg_start=1;
end
if bkg(2) > max(size(spectrum))
     bkg(2)=max(size(spectrum));
end
varargout{1}=map2colvec(bkg);
if nargout>1
    % Set indices for zooming in on spectrum
    title('Select region of interest')
    temp=ginput(2);
    crop(1)=uint16(round(min(temp(1,1),temp(2,1))));
    crop(2)=uint16(round(max(temp(1,1),temp(2,1))));
    if crop(1) < 1
        crop(1)=1;
    end
    if crop(2) > max(size(spectrum))
        crop(2)=max(size(spectrum));
    end
    plot((1:length(spectrum(crop(1):crop(2))))+double(crop(1))-1,spectrum(crop(1):crop(2)),'k')
    varargout{2}=map2colvec(crop);
end
%     spectrum=spectrum(crop(1):crop(2));
if nargout==3
% Set indices for boxcar
    title('Choose start and end of boxcar')
    temp=ginput(2);
    boxcar(1)=uint16(round(min(temp(1,1),temp(2,1))));
    boxcar(2)=uint16(round(max(temp(1,1),temp(2,1))));
    if boxcar(1) < 1
        boxcar(1)=1;
    end
    if boxcar(2) > max(size(spectrum))
        boxcar(2)=max(size(spectrum));
    end
    varargout{3}=map2colvec(boxcar);
end
close(hndl1)
% varargout{3}=map2colvec((spectrum(crop(1):crop(2))-mean(spectrum(crop(1):crop(2))))/max(spectrum(crop(1):crop(2))-mean(spectrum(crop(1):crop(2)))));
end
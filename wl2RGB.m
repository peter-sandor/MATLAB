function RGBout = wl2RGB(wl)

% Wavelength to RGB color conversion from
% http://www.physics.sfasu.edu/astro/color/spectra.html
% wl in nm
for ind1=1:length(wl)
    if ((wl(ind1)>=380.)&&(wl(ind1)<440.))
        RGBout(ind1,:) = [-1.*(wl(ind1)-440)/(440-380) 0 1];
        if wl(ind1)<=420
            RGBout(ind1,:)=(0.3+0.7*(wl(ind1)-380)/(420-380))*RGBout(ind1,:); 
        end
    elseif ((wl(ind1)>=440.)&&(wl(ind1)<490.))
        RGBout(ind1,:) = [0 (wl(ind1)-440.)/(490.-440.) 1];
    elseif ((wl(ind1)>=490.)&&(wl(ind1)<510.))
        RGBout(ind1,:) = [0 1 -1.*(wl(ind1)-510.)/(510.-490.)];
    elseif ((wl(ind1)>=510.)&&(wl(ind1)<580.))
        RGBout(ind1,:) = [(wl(ind1)-510.)/(580.-510.) 1 0];
    elseif ((wl(ind1)>=580.)&&(wl(ind1)<645.))
        RGBout(ind1,:) = [1 -1.*(wl(ind1)-645.)/(645.-580.) 0];
    elseif ((wl(ind1)>=645)&&(wl(ind1)<=780.))
        RGBout(ind1,:) = [1 0 0];
        if wl(ind1)>=700
            RGBout(ind1,:)=(0.3+0.7*(780-wl(ind1))/(780-700))*RGBout(ind1,:); 
        end
    else
        RGBout(ind1,:)=[0 0 0];
    end
end

end
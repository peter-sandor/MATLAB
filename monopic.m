function pic = monopic(pic0)
% converts a given 8bit picture to a monochromatic single precison pciture

mono = round((3 - length(size(pic0)))/2);
if mono
        pic = double(pic0);
    else
        pic = double(pic0(:,:,1));    
end;




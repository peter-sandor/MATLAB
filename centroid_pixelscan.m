function varargout = centroid_pixelscan(Image, Threshold)
%check if an image has exactly one hit and register its coord

%check each pixel against the 'Threshold'. If above threshold, check
%adjacent 4 pixels. If two of them, upper left or lower right, then
%register the its coord.

%tic
PixeList = [];
for j = 2:359
    for k = 2:359
        if Image(j,k)>Threshold
            if Image(j-1,k)>Threshold && Image(j,k-1)>Threshold
                PixeList = cat(1,PixeList, [j k]);
            elseif Image(j+1,k)>Threshold && Image(j,k+1)>Threshold
                PixeList = cat(1,PixeList, [j k]);
            end
        end
    end
end


%calculate the center of mass of all registered pixels and output it.

MomentJ =0;
MomentK =0;
Mass = 0;
if length(PixeList)>4 && (PixeList(end,1)-PixeList(1,1)+abs(PixeList(1,2)-PixeList(end,2)))<16
    for m = 1:length(PixeList)
        MomentJ = MomentJ + Image(PixeList(m,1),PixeList(m,2))*PixeList(m,1);
        MomentK = MomentK + Image(PixeList(m,1),PixeList(m,2))*PixeList(m,2);
        Mass = Mass + Image(PixeList(m,1),PixeList(m,2));
    end
    Coord = [MomentK/Mass MomentJ/Mass];
    Hit=1;
else
    Hit=0;
    Coord=[-1 -1];
end
varargout{1}=Hit;
varargout{2}=Coord;
%t=toc
end


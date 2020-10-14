function coords = calculate_coordinates(S,minlistlength)

hits=0;
coords=[];
for ind1=1:length(S)
    if length(S{ind1})>minlistlength % checks whether each list contains a number of pixels that is reasonable large
                                     % this helps filter out single pixel "defects" with large illumination value on them, which may result in false hits
        hits=hits+1;
        coords(hits,:)=[0 0];
        for ind2=1:length(S{ind1})
            coords(hits,:)=coords(hits,:)+S{ind1}(ind2,3)*S{ind1}(ind2,1:2);
        end
        coords(hits,:)=coords(hits,:)/sum(S{ind1}(:,3));
    end
end
end
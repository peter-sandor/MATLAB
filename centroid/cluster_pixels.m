function S = cluster_pixels(pts,radius,minlistlength)

indices=1;
for ind1=1:size(pts,1)-1
    if pts(ind1+1,2)-pts(ind1,2)>=radius
        indices=[indices; ind1+1];
    end
end

ind2=1;
for ind1=1:length(indices)-1
    temp=pts(indices(ind1):indices(ind1+1)-1,:);
    if max(temp(:,1))-min(temp(:,1))>5*radius
        temp2=sortrows(temp);
        indices2=1;
        for ind3=1:size(temp2,1)-1
            if temp2(ind3+1,1)-temp2(ind3,1)>=radius
                indices2=[indices2; ind3+1];
            end
        end
        for ind3=1:length(indices2)-1
            if (indices2(ind3+1)-indices2(ind3))>minlistlength
                S{ind2}=temp2(indices2(ind3):indices2(ind3+1)-1,:);
                ind2=ind2+1;
            end
        end
        if (size(temp2,1)-indices2(end))>minlistlength
            S{ind2}=temp2(indices2(end):end,:);
            ind2=ind2+1;
        end
    else
        if (indices(ind1+1)-indices(ind1))>minlistlength
            S{ind2}=pts(indices(ind1):indices(ind1+1)-1,:);
            ind2=ind2+1;
        end
    end
end

temp=pts(indices(end):end,:);
if max(temp(:,1))-min(temp(:,1))>5*radius
    temp2=sortrows(temp);
    indices2=1;
    for ind3=1:size(temp2,1)-1
        if temp2(ind3+1,1)-temp2(ind3,1)>=radius
            indices2=[indices2; ind3+1];
        end
    end
    for ind3=1:length(indices2)-1
        if (indices2(ind3+1)-indices2(ind3))>minlistlength
            S{ind2}=temp2(indices2(ind3):indices2(ind3+1)-1,:);
            ind2=ind2+1;
        end
    end
    if (size(temp2,1)-indices2(end))>minlistlength
        S{ind2}=temp2(indices2(end):end,:);
        ind2=ind2+1;
    end
else
    if (size(pts,1)-indices(end))>minlistlength
        S{ind2}=pts(indices(end):end,:);
        ind2=ind2+1;
    end
end

end
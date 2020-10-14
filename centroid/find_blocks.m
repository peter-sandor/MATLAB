for ind1=1:size(Mconnect,1)-1
%     temp(ind1)=sum(double(Mconnect(ind1:end,ind1)),1)+sum(double(Mconnect(ind1+1,1:ind1)),2)-1;
%     temp2(ind1)=-1;
%     if ind1<=size(Mconnect,1)/2
%         bnd=ind1;
%     else
%         bnd=size(Mconnect,1)-ind1;
%     end
%     for ind2=1:min(100,bnd)
%         temp2(ind1)=temp2(ind1)+double(Mconnect(ind1-ind2+1,ind1+ind2-1));
%     end
    if ind1<=size(Mconnect,1)/2
        range=min(ind1-1,50);
    else
        range=min(size(Mconnect,1)-ind1,50);
    end
    temp5(ind1)=sum(sum(double(Mconnect(ind1-range:ind1,ind1+1:ind1+range))));
%     disp([num2str(ind1) ' ' num2str(range)]);
end

sumlogic=true([length(summed) 1]);
sumlogic(summed<=5)=false;

flag=false;
ind2=1;
for ind1=1:length(sumlogic)
    if sumlogic(ind1) && ~flag % start of range of pixels for a hit
        flag=true;
        indices(ind2,1)=ind1;
        
    elseif ~sumlogic(ind1) && flag % end of range of pixels
        indices(ind2,2)=ind1;
        flag=false;
        ind2=ind2+1;
    end
end

for ind1=1:size(indices,1)
    S{ind1}=pts(indices(ind1,1):indices(ind1,2),:);
end
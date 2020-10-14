function indices = presort_pixels(Mconnect,filtersize)

summed=zeros([size(Mconnect,1)-1 1]);
for ind1=1:size(Mconnect,1)-1
    if ind1<=size(Mconnect,1)/2
        range=min(ind1-1,filtersize);
    else
        range=min(size(Mconnect,1)-ind1,filtersize);
    end
    summed(ind1)=sum(sum(double(Mconnect(ind1-range:ind1,ind1+1:ind1+range))));
%     disp([num2str(ind1) ' ' num2str(range)]);
end

sumlogic=true([length(summed) 1]);
sumlogic(summed<=5)=false;

% indices=zeros([length(sumlogic) 2]);
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
if ~exist('indices','var');
    indices=[];
end
end
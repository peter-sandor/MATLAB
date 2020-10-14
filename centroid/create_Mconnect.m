function Mconnect = create_Mconnect(pts,radius,filtersize)
% Create connection matrix with size M x M, logical type. Two pixels are connected
% if their distance is smaller than the value of 'radius'.
M=size(pts,1);
% ind2=1;
% ind_keep=[];
Mconnect=zeros([M M]);
for ind1=1:M
    temp=((pts(ind1,1)-pts(:,1)).^2+(pts(ind1,2)-pts(:,2)).^2)<=radius^2;
%     if sum(temp)>minlistlength % filter elements of the list additionally for standalone pixels with no neighbours that survived the thresholding. These are surely not valid hits.
        Mconnect(:,ind1)=logical(temp);
%         ind2=ind2+1;
%         ind_keep=[ind_keep; ind1];
%     end
end
%         Mconnect=Mconnect(ind_keep,:);
%         pts=pts(ind_keep,:);
end
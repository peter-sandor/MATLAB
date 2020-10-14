function S = cluster_pixels_M(pts,Mconnect)

ind1=1;
while ~isempty(Mconnect)
    % find first-order connections
    member0=logical(Mconnect(:,1)); % index vector (column)
    % find second-order connections
    member1=any(Mconnect(member0,:),1);
    member1=map2colvec(member1); % index vector (column)
    member_diff=xor(member0,member1);
    Mconnect(member0,:)=false;
    Mconnect(:,map2rowvec(member0))=false;
    member_tot=member1;
    % check and follow up on higher-order connections
    while any(member_diff)
        member0=member1;
        member1=any(Mconnect(member_diff,:),1);
        member1=map2colvec(logical(member1)); % index vector (column)
        member_diff=xor(member0,member1);
        member_tot=or(member_tot,member1);
        Mconnect(member0,:)=false;
        Mconnect(:,map2rowvec(member0))=false;
    end
    S{ind1}=pts(member_tot,:); % 'member_tot' contains the indices for all the connected pixels, here create a separate list for them. 
    pts(member_tot,:)=[]; % Delete elements that have already been taken into account. This speeds up the sorting.
    Mconnect(member_tot,:)=[];
    Mconnect(:,map2rowvec(member_tot))=[];
    ind1=ind1+1;
end

end
function out = filter_to_range(in,lim)
% This function deletes the ordered M-tuples in the input array 'in' which have at least element beyond a specified range.
% 'in' has a size of NxM, so the M-tuples are row vectors.
% Ranges for the different elements in the tuples can be different, and are
% contained in the columns of 'lim', such that:
% lim(:,1) --> in(:,1) and lim(:,2) --> in(:,2), etc.

N=size(in,1);
if ~isempty(lim)
    ind_vec=ones([N 1]);
    for ind1=1:size(in,2)
        ind_vec=ind_vec.*map2colvec(in(:,ind1)<=max(lim(:,ind1))).*map2colvec(in(:,ind1)>=min(lim(:,ind1)));
    end
    out=in(logical(ind_vec),:);
else
    out=in;
end
end
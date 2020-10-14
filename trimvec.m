function varargout = trimvec(totrim,reference)
% Trim a vector to make sure it doesn't have elements smaller or bigger
% than that of a reference vector. Useful for consolidating delay axes.
index_start=max(vec2ind(totrim<min(reference)))+1;
index_end=min(vec2ind(totrim>max(reference)))-1;
if isempty(index_start)
    index_start=1;
end
if isempty(index_end)
    index_end=length(totrim);
end
varargout{1}=totrim(index_start:index_end);
varargout{2}=[index_start index_end];
end
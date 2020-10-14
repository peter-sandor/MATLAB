function vec_out = find4thcorner(vec_in)
% Given 3 vectors pointing to the corners of a rectangle, find the 4th one,
% and insert it to the array that contains the others. Preserve the order
% of the vectors in the array. (The order defines the orientation of the
% surface normal!)
for ind1 = 1:3
    vec_shift = circshift(vec_in,[ind1-1 0]);
    dot_product = map2rowvec(vec_shift(2,:)-vec_shift(1,:))*map2colvec(vec_shift(3,:)-vec_shift(1,:));
    if dot_product == 0
        vec_out = circshift(vec_shift,[1 0]);
        vec_out(4,:) = vec_shift(2,:) + vec_shift(3,:) - vec_shift(1,:);
    end
end
end
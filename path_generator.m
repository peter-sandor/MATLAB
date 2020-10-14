function path_out = path_generator(positions,start_pos)

for ind0=1:length(positions)
    N(ind0) = length(positions{ind0});
end
temp = cumprod(N);
N_points = temp(end);
path_out = zeros([N_points, length(N)]);
% path_out(1,:) = start_pos;

for ind0 = 1:N_points
    ind2 = ceil(ind0/N(1));
    ind3 = ceil(ind0/N(1)/N(2));
    if mod(ind2,2) == 0
        ind1 = N(1) - mod(ind0-1,N(1));
    else
        ind1 = mod(ind0-1,N(1)) + 1;
    end
    path_out(ind0,:) = [positions{1}(ind1) positions{2}(ind2) positions{3}(ind3)];
end
end
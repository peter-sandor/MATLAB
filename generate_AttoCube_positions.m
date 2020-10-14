function varargout = generate_AttoCube_positions(params,filename)
% This function generates the position file for scans with the AttoCube
% params: 3x3 array; k-th row is of the form: [start_position(k) end_position(k)
% stepsize(k)] in micrometers for the axis 'k' (A1--> k=1, etc.)

for ind1 = 1:3
    positions{ind1} = [params(ind1,1):params(ind1,3):params(ind1,2)]*1000;
end

% path_out = combine(positions);
start_pos = [positions{1}(1) positions{2}(1) positions{3}(1)];
path_out = path_generator(positions,start_pos);
varargout{1} = path_out;
dlmwrite(filename,path_out,'delimiter','\t','precision','%.0f');
end
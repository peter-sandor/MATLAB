function varargout = construct_image_from_scan(filename)

data_dot = load(filename);
data1 = data_dot;
data1(:,1:3) = roundP(data_dot(:,1:3),0);
N = size(data1,1);

% fix irregularities in stepsize in A1
pos_A1_round = roundP(data_dot(:,1),0);
% vec_stepsize_A1 = 0:max(diff(unique(pos_A1_round)));
vec_stepsize_A1 = 0:max(abs(diff(pos_A1_round)));
hist1 = histcounts(abs(diff(pos_A1_round)),vec_stepsize_A1+0.5);
stepsize_A1 = vec_stepsize_A1(vec2ind(hist1 == max(hist1))+1);
if mod(max(pos_A1_round)-min(pos_A1_round),stepsize_A1)==0
    Npos(1) = (max(pos_A1_round)-min(pos_A1_round))/stepsize_A1+1;
    pos_vec{1} = min(pos_A1_round):stepsize_A1:max(pos_A1_round);
    for ind1 = 1:size(data1,1)
        ind_insert = vec2ind(abs(data1(ind1,1)-pos_vec{1})<=stepsize_A1/2);
        data1(ind1,1) = pos_vec{1}(ind_insert);
    end
else
    Npos(1) = 0;
end

% fix irregularities in stepsize in A2
pos_A2_round = roundP(data_dot(:,2),0);
% vec_stepsize_A2 = 0:max(diff(unique(pos_A2_round)));
vec_stepsize_A2 = 0:max(abs(diff(pos_A2_round)));
hist2 = histcounts(abs(diff(pos_A2_round)),vec_stepsize_A2+0.5);
stepsize_A2 = vec_stepsize_A2(vec2ind(hist2 == max(hist2))+1);
if mod(max(pos_A2_round)-min(pos_A2_round),stepsize_A2)==0
    Npos(2) = (max(pos_A2_round)-min(pos_A2_round))/stepsize_A2+1;
    pos_vec{2} = min(pos_A2_round):stepsize_A2:max(pos_A2_round);
    for ind1 = 1:size(data1,1)
        ind_insert = vec2ind(abs(data1(ind1,2)-pos_vec{2})<=stepsize_A2/2);
        data1(ind1,2) = pos_vec{2}(ind_insert);
    end
else
    Npos(2) = 0;
end

% pos_vec{1} = unique(data1(:,1));
% pos_vec{2} = unique(data1(:,2));
% Npos(2) = length(pos_vec{1});
% Npos(1) = length(pos_vec{2});

pic1 = zeros([Npos(2) Npos(1)]);
for ind1 = 1:N
    ind_pos1 = vec2ind(pos_vec{1} == data1(ind1,1));
    ind_pos2 = vec2ind(pos_vec{2} == data1(ind1,2));
    pic1(ind_pos2,ind_pos1) = data1(ind1,4);
end
pic2 = pic1 - min(min(pic1));
varargout{1} = pic2;
varargout{2} = pos_vec;
end
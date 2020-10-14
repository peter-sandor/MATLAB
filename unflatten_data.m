function varargout = unflatten_data(data1)
% data1 is an NxM array
% This code is not operational!
N = size(data1);

for ind1 = 1:N(2)-1
    pos_vec{ind1} = unique(data1(:,ind1));
    Npos(ind1) = length(pos_vec{ind1});
end

% pic1 = zeros(Npos);
pic1 = reshape(data1(:,end),Npos);

% for ind2 = 1:N(1)
%     for ind3 = 1:N(2)-1
%         ind_pos(ind3) = vec2ind(pos_vec{ind3} == data1(ind2,ind3));
%         pic1(ind_pos2,ind_pos1) = data1(ind2,end-1);
%     end
% end
varargout{1} = pic1;
varargout{2} = pos_vec;
end
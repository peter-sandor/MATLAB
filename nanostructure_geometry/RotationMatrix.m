function R = RotationMatrix(vec_dir)
% This function calculates the rotation matrix for a rotation that takes
% the unit vector in the 'z' direction ([0 0 1]) and rotates it to point in the
% direction of input vector 'vec_dir', for which the components are [x y
% z].
% From section 9.2 of Modelling CPV by Ian Richard Cole (https://dspace.lboro.ac.uk/2134/18050)

% find axis vector 'vec_axis' (it is perpendicular to 'vec_dir' and z-axis)
vec_axis = cross([0 0 1],vec_dir);
if norm(vec_axis)~=0
    vec_axis = vec_axis/norm(cross([0 0 1],vec_dir));
end
theta = acos(map2rowvec([0 0 1])*map2colvec(vec_dir)); % rotation angle
R(1,1) = cos(theta) + vec_axis(1)^2*(1-cos(theta));
R(1,2) = vec_axis(1)*vec_axis(2)*(1-cos(theta)) - vec_axis(3)*sin(theta);
R(1,3) = vec_axis(1)*vec_axis(3)*(1-cos(theta)) + vec_axis(2)*sin(theta);
R(2,1) = vec_axis(2)*vec_axis(1)*(1-cos(theta)) + vec_axis(3)*sin(theta);
R(2,2) = cos(theta) + vec_axis(2)^2*(1-cos(theta));
R(2,3) = vec_axis(2)*vec_axis(3)*(1-cos(theta)) - vec_axis(1)*sin(theta);
R(3,1) = vec_axis(3)*vec_axis(1)*(1-cos(theta)) - vec_axis(2)*sin(theta);
R(3,2) = vec_axis(3)*vec_axis(2)*(1-cos(theta)) + vec_axis(1)*sin(theta);
R(3,3) = cos(theta) + vec_axis(3)^2*(1-cos(theta));
end
function [point_prj,surf_norm] = nanostructure_analysis(points_in)
% Define nanostructure
side_x = 235; % [nm]
side_y = 135;
side_z = 40;
radius = 15;

element(1).type = 'rectangle';
element(1).vec = [-(side_x/2-radius) -(side_y/2) side_z-radius;
                  -(side_x/2-radius) -(side_y/2) 0;
                  +side_x/2-radius -(side_y/2) 0;];
           
element(2).type = 'rectangle';
element(2).vec = [+side_x/2 -(side_y/2-radius) side_z-radius;
               +side_x/2 -(side_y/2-radius) 0;
               +side_x/2 +side_y/2-radius 0];

element(3).type = 'rectangle';
element(3).vec = [+side_x/2-radius +side_y/2 side_z-radius;
               +side_x/2-radius +side_y/2 0;
               -(side_x/2-radius) +side_y/2 0];

element(4).type = 'rectangle';
element(4).vec = [-(side_x/2) +side_y/2-radius side_z-radius;
               -(side_x/2) +side_y/2-radius 0;
               -(side_x/2) -(side_y/2-radius) 0];

element(5).type = 'rectangle';
element(5).vec = [-(side_x/2-radius) -(side_y/2-radius) side_z;
               +side_x/2-radius -(side_y/2-radius) side_z;
               +side_x/2-radius +side_y/2-radius side_z;];

element(6).type = 'sphere';
element(6).vec = [-(side_x/2-radius) -(side_y/2-radius) side_z-radius];
element(6).R = radius;

element(7).type = 'sphere';
element(7).vec = [+side_x/2-radius -(side_y/2-radius) side_z-radius];
element(7).R = radius;

element(8).type = 'sphere';
element(8).vec = [+side_x/2-radius +side_y/2-radius side_z-radius];
element(8).R = radius;

element(9).type = 'sphere';
element(9).vec = [-(side_x/2-radius) +side_y/2-radius side_z-radius];
element(9).R = radius;

element(10).type = 'cylinder';
element(10).vec = [-(side_x/2-radius) -(side_y/2-radius) side_z-radius;
                +side_x/2-radius -(side_y/2-radius) side_z-radius];
element(10).R = radius;

element(11).type = 'cylinder';
element(11).vec = [+side_x/2-radius -(side_y/2-radius) side_z-radius;
                +side_x/2-radius +side_y/2-radius side_z-radius];
element(11).R = radius;

element(12).type = 'cylinder';
element(12).vec = [+side_x/2-radius +side_y/2-radius side_z-radius;
                -(side_x/2-radius) +side_y/2-radius side_z-radius];
element(12).R = radius;

element(13).type = 'cylinder';
element(13).vec = [-(side_x/2-radius) +side_y/2-radius side_z-radius;
                -(side_x/2-radius) -(side_y/2-radius) side_z-radius];
element(13).R = radius;

element(14).type = 'cylinder';
element(14).vec = [-(side_x/2-radius) -(side_y/2-radius) 0;
                -(side_x/2-radius) -(side_y/2-radius) side_z-radius;];
element(14).R = radius;

element(15).type = 'cylinder';
element(15).vec = [+(side_x/2-radius) -(side_y/2-radius) 0;
                +(side_x/2-radius) -(side_y/2-radius) side_z-radius;];
element(15).R = radius;

element(16).type = 'cylinder';
element(16).vec = [+(side_x/2-radius) +(side_y/2-radius) 0;
                +(side_x/2-radius) +(side_y/2-radius) side_z-radius;];
element(16).R = radius;

element(17).type = 'cylinder';
element(17).vec = [-(side_x/2-radius) +(side_y/2-radius) 0;
                -(side_x/2-radius) +(side_y/2-radius) side_z-radius];
element(17).R = radius;
%% Find points and surface normals on the nanostructure
for ind1 = 1:size(points_in,1)
    norms = [];
    for ind2 = 1:length(element)
        [temp1 temp2] = project2surface(points_in(ind1,:),element(ind2));
        if ~isempty(temp1)
            norms(ind2) = norm(points_in(ind1,:)-temp1);
            vec_proj(ind2,:) = temp1;
            vec_norm(ind2,:) = temp2;
        else
            norms(ind2) = Inf;
            vec_proj(ind2,:) = [0 0 0];
            vec_norm(ind2,:) = [0 0 0];
        end
    end
    [a,ind_norm] = min(norms);
    point_prj(ind1,:) = vec_proj(ind_norm,:);
    surf_norm(ind1,:) = vec_norm(ind_norm,:);
end
%% Plot
figure;hax = axes;
for ind1 = 1:length(element)
    draw_element(hax,element(ind1))
end
axis equal;
hold(hax,'on');
for ind1 = 1:size(points_in,1)
%     quiver3(hax,point_prj(ind1,1),point_prj(ind1,2),point_prj(ind1,3),10*surf_norm(ind1,1),10*surf_norm(ind1,2),10*surf_norm(ind1,3),'Color','red','Linewidth',3)
    quiver3(hax,points_in(ind1,1),points_in(ind1,2),points_in(ind1,3),10*surf_norm(ind1,1),10*surf_norm(ind1,2),10*surf_norm(ind1,3),'Color','blue','Linewidth',3)
end
hold(hax,'off');
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')
end
function draw_element(hax,element)
Ndraw  = 25;
clr = [0.8549,0.6471,0.1255];
hold(hax,'on');
if strcmp(element.type,'rectangle')
    vecs = find4thcorner(element.vec);
    patch(hax,'Xdata',vecs(:,1),'Ydata',vecs(:,2),'Zdata',vecs(:,3),'FaceColor',clr,'EdgeColor','black')
elseif strcmp(element.type,'sphere')
    [X,Y,Z] = sphere(Ndraw);
    surf(hax,'Xdata',X*element.R + element.vec(1),'Ydata',Y*element.R + element.vec(2),'Zdata',Z*element.R + element.vec(3),'FaceColor',clr,'EdgeColor','black')
elseif strcmp(element.type,'cylinder')
    [X,Y,Z] = cylinder(ones([Ndraw 1]),Ndraw-1);
    vec_dir = (element.vec(2,:)-element.vec(1,:))/norm(element.vec(2,:)-element.vec(1,:));
    X2 = element.R * X;
    Y2 = element.R * Y;
    Z2 = norm(element.vec(2,:)-element.vec(1,:)) * Z;
    R_rot = RotationMatrix(vec_dir);
    for ind1 = 1:Ndraw
        for ind2 = 1:Ndraw
            temp = R_rot * [X2(ind1,ind2); Y2(ind1,ind2); Z2(ind1,ind2)];
            X3(ind1,ind2) = temp(1) + element.vec(1,1);
            Y3(ind1,ind2) = temp(2) + element.vec(1,2);
            Z3(ind1,ind2) = temp(3) + element.vec(1,3);
        end
    end
    surf(hax,'Xdata',X3,'Ydata',Y3,'Zdata',Z3,'FaceColor',clr,'EdgeColor','black');
end
hold(hax,'off');
end
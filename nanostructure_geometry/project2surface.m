function [vec_proj surf_norm] = project2surface(vec_in,element)

if strcmp(element.type,'rectangle')
    vecs = find4thcorner(element.vec);
    surf_norm = cross(vecs(2,:)-vecs(1,:),vecs(3,:)-vecs(1,:))/norm(cross(vecs(2,:)-vecs(1,:),vecs(3,:)-vecs(1,:)));
    t = map2rowvec(vecs(1,:)) * map2colvec(surf_norm);
    vec_proj = vec_in - (map2rowvec(vec_in - t*surf_norm)*map2colvec(surf_norm))*surf_norm;
    a = map2rowvec(vec_proj - vecs(1,:))*map2colvec(vecs(2,:)-vecs(1,:))/norm(vecs(2,:)-vecs(1,:))^2;
    b = map2rowvec(vec_proj - vecs(1,:))*map2colvec(vecs(4,:)-vecs(1,:))/norm(vecs(4,:)-vecs(1,:))^2;
    if a<0 || a>1 || b<0 || b>1
        vec_proj = [];
    end
elseif strcmp(element.type,'sphere')
    surf_norm = (vec_in - element.vec)/norm(vec_in - element.vec);
    vec_proj = element.vec + surf_norm*element.R;
elseif strcmp(element.type,'cylinder')
    s = map2rowvec(vec_in - element.vec(1,:))*map2colvec(element.vec(2,:)-element.vec(1,:))/norm(element.vec(2,:)-element.vec(1,:))^2;
    if s<0 && norm(vec_in-(s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:))) < element.R
        surf_norm = (element.vec(1,:)-element.vec(2,:))/norm(element.vec(1,:)-element.vec(2,:));
        t = map2rowvec(element.vec(1,:)) * map2colvec(surf_norm);
        vec_proj = vec_in - (map2rowvec(vec_in - t*surf_norm)*map2colvec(surf_norm))*surf_norm;
    elseif s>1 && norm(vec_in-(s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:))) < element.R
        surf_norm = (element.vec(2,:)-element.vec(1,:))/norm(element.vec(2,:)-element.vec(1,:));
        t = map2rowvec(element.vec(2,:)) * map2colvec(surf_norm);
        vec_proj = vec_in - (map2rowvec(vec_in - t*surf_norm)*map2colvec(surf_norm))*surf_norm;
    elseif s>=0 && s<=1 
        surf_norm = (vec_in - (s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:)))/norm(vec_in - (s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:)));
        vec_proj = element.vec(1,:) + s*(element.vec(2,:)-element.vec(1,:)) + element.R*(vec_in - (s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:)))/norm(vec_in - (s*(element.vec(2,:)-element.vec(1,:)) + element.vec(1,:)));
    else
        surf_norm = [];
        vec_proj = [];
    end
end

end
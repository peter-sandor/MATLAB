function line_refl = calc_reflection(line_in,object_in,point_in)

PREC=1e-10;
x=point_in(1);
y=point_in(2);
if  isfield(object_in,'A')
    line_surf=object_in;
    % calculate point of reflection
    surf_normal=[-line_surf.k(2); line_surf.k(1)];
elseif isfield(object_in,'R')
	surf_normal=([object_in.center(1)-x; object_in.center(2)-y])/sqrt((object_in.center(1)-x)^2+(object_in.center(2)-y)^2);
end
if (line_in.k(1)*surf_normal(1)+line_in.k(2)*surf_normal(2))>0
    surf_normal=-surf_normal;
end
temp=line_in.k+2*abs(map2rowvec(line_in.k)*map2colvec(surf_normal))*surf_normal;
if abs(temp(2))<=PREC
    temp(2)=0;
end
if abs(temp(1))<=PREC
    temp(1)=0;
    line_refl.A=Inf;
    line_refl.B=x;
else
    line_refl.A=temp(2)/temp(1);
    line_refl.B=y-line_refl.A*x;
end
line_refl.k=temp;
line_refl.pstart=point_in;
line_refl.pend=[];
end
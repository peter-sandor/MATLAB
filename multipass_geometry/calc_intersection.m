function [points_out]=calc_intersection(object1,object2)
% straight = [A B] for y=A*x+B; if A=inf, then line_str(2)=x coord.
% curved = [x0 y0 R] for R^2=(x-x0)^2+(y-y0)^2
if isfield(object1,'A')
    ID1='line';
    A1=object1.A;
    B1=object1.B;
elseif isfield(object1,'R')
    ID1='circle';
    x0=object1.center(1);
    y0=object1.center(2);
    R=object1.R;
end
if isfield(object2,'R')
    ID2='circle';
    x0=object2.center(1);
    y0=object2.center(2);
    R=object2.R;
elseif isfield(object2,'A')
    ID2='line';
    if strcmp(ID1,'circle')
        A1=object2.A;
        B1=object2.B;
    else
        A2=object2.A;
        B2=object2.B;
    end
end

if strcmp(ID1,'circle') || strcmp(ID2,'circle') % calculate intersection between line and circle
    if A1==0 % line is parallel to x-axis
        A_par=1;
        B_par=-2*x0;
        C_par=x0^2-R^2+(B1-y0)^2;
        if B_par^2>4*A_par*C_par
            x=[(-B_par+sqrt(B_par^2-4*A_par*C_par))/2/A_par; (-B_par-sqrt(B_par^2-4*A_par*C_par))/2/A_par];
            y=[B1; B1];
        elseif B_par^2==4*A_par*C_par
            y=B1;
            x=-B_par/2/A_par;
        else
            y=[];
            x=[];
        end
    elseif isinf(A1) % line is parallel to y-axis
        A_par=1;
        B_par=-2*y0;
        C_par=y0^2-R^2+(B1-x0)^2;
        if B_par^2>4*A_par*C_par
            y=[(-B_par+sqrt(B_par^2-4*A_par*C_par))/2/A_par; (-B_par-sqrt(B_par^2-4*A_par*C_par))/2/A_par];
            x=[B1; B1];
        elseif B_par^2==4*A_par*C_par
            x=B1;
            y=-B_par/2/A_par;
        else
            y=[];
            x=[];
        end
    else % general case, line is not parallel to either axes
        A_par=1+1/A1^2;
        B_par=-2*(y0+B1/A1^2+x0/A1);
        C_par=y0^2+(B1/A1)^2-R^2+2*x0*B1/A1+x0^2;
        if B_par^2>4*A_par*C_par
            y=[(-B_par+sqrt(B_par^2-4*A_par*C_par))/2/A_par; (-B_par-sqrt(B_par^2-4*A_par*C_par))/2/A_par];
            x=(y-B1)/A1;
        elseif B_par^2==4*A_par*C_par
            y=-B_par/2/A_par;
            x=(y-B1)/A1;
        else
            y=[];
            x=[];
        end
    end
elseif strcmp(ID1,'line') && strcmp(ID2,'line') % calculate intersection between two straight lines
    if A1==A2 % two straight lines are parallel
        x=[];
        y=[];
    else
        if isinf(A1) % one of the lines is parallel to y-axis
            x=B1;
            y=A2*x+B2;
        elseif isinf(A2) % other line is parallel to y-axis
            x=B2;
            y=A1*x+B1;
        else % general case, neither of the lines is parallel to y-axis and not parallel to each other
            x=(B1-B2)/(A2-A1);
            y=A1*x+B1;
        end
    end
end
points_out=[map2colvec(x) map2colvec(y)];
end
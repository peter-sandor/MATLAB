function varargout=draw_braces(dims,orient,dratio)

% dims = [x0 y0 xsize ysize]
x0=dims(1);
y0=dims(2);
xsize=dims(3);
ysize=dims(4);
Nb=100;

if strcmp(orient,'right')
    pts1=[x0+xsize y0;
        x0 y0;
        x0+xsize y0+ysize*dratio;
        x0 y0+ysize*dratio];
    pts2=[x0+xsize y0+ysize;
        x0 y0+ysize;
        x0+xsize y0+ysize*dratio;
        x0 y0+ysize*dratio];
elseif strcmp(orient,'left')
    pts1=[x0 y0;
        x0+xsize y0;
        x0 y0+ysize*dratio;
        x0+xsize y0+ysize*dratio];
    pts2=[x0 y0+ysize;
        x0+xsize y0+ysize;
        x0 y0+ysize*dratio;
        x0+xsize y0+ysize*dratio];
end

curve1=Bezier_curve(pts1,Nb);
curve2=Bezier_curve(pts2,Nb);
hndl(1)=line(curve1(:,1),curve1(:,2));
hndl(2)=line(curve2(:,1),curve2(:,2));
varargout{1}=hndl;
end
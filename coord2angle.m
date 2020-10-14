function out = coord2angle(in)

% in this code, alpha=0 corresponds to the positive y-axis (associated with
% row index), and increasing alpha takes us around counterclockwise.
x=in(1);
y=in(2);
if y<0
    out = atan(-x/y)+pi;
elseif y>0
    if x<=0
        out = atan(-x/y);
    elseif x>0
        out = atan(-x/y)+2*pi;
    end
elseif y==0
    if x<0
        out = pi/2;
    elseif x>0
        out = 3*pi/2;
    else 
        out=[];
    end
end

end
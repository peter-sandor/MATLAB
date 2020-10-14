function out = coord2angle_original(in)

% in this code, alpha=0 corresponds to positive x-axis (associated with
% column index), and increasing alpha takes us around counterclockwise.
x=in(1);
y=in(2);
if x<0
    out = atan(y/x)+pi;
elseif x>0
    if y>=0
        out = atan(y/x);
    elseif y<0
        out = atan(y/x)+2*pi;
    end
elseif x==0
    if y>0
        out = pi/2;
    elseif y<0
        out = 3*pi/2;
    else 
        out=[];
    end
end

end
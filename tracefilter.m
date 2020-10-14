function out=tracefilter(filename,thrs,ROI,varargin)
trace=load(filename);
if nargin==4;
    invert=varargin{1};
else 
    invert=0;
end
if invert==0
    if max(trace(ROI(1):ROI(2)))>=thrs
        out=1;
    else
        out=0;
    end
elseif invert==1
    if min(trace(ROI(1):ROI(2)))<=thrs
        out=1;
    else
        out=0;
    end
end
end
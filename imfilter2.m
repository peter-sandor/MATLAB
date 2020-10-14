function out=imfilter2(filename,thrs,varargin)
pic=loadbinaryimage2(filename);
if nargin==3;
    region=varargin{1};
else region=[1 size(pic,1);1 size(pic,2)];
end
if max(max(pic(region(1,1):region(1,2),region(2,1):region(2,2)))) >= thrs
    out=1;
else out=0;
end
end
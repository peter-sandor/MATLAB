function compactify(varargin)
if nargin == 1
    hfig = varargin{1};
else
    hfig = gcf;
end
set(hfig,'outerposition',[0.25 0.15 0.5 0.65])
end


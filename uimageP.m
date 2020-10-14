function varargout = uimage(varargin)
%UIMAGEP  Display image with uneven axis.
%   UIMAGE(X,Y,C) displays matrix C as an image, using the vectors X and
%   Y to specify the X and Y coordinates. X and Y may be unevenly spaced
%   vectors, but must be increasing. The size of C must be LENGTH(Y)*
%   LENGTH(X). (Most probably you'll want to display C' instead of C).
%
%   Contrary to Matlab's original IMAGE function, here the vectors X and Y
%   do not need to be linearly spaced. Whereas IMAGE linearly interpolates
%   the X-axis between X(1) and X(end), ignoring all other values (idem
%   for Y), UIMAGE allows for X and/or Y to be unevenly spaced vectors, by

error(nargchk(3,inf,nargin));

if nargin==1
    z=varargin;
    x=1:size(z,2);
    y=1:size(z,1);
elseif nargin==3
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
end

% calculate vertices
vx = [(x(1)-x(2))/2+x(1) (x(1:end-1) + diff(x)/2) (x(end)-x(end-1))/2+x(end)];
vy = [(y(1)-y(2))/2+y(1) (y(1:end-1) + diff(y)/2) (y(end)-y(end-1))/2+y(end)];

hax0=axes;
for ind1=1:length(vx)-1
    for ind2=1:length(vy)-1
        hpat(ind2,ind1) = patch([vx(ind1) vx(ind1+1) vx(ind1+1) vx(ind1)],[vy(ind2) vy(ind2) vy(ind2+1) vy(ind2+1)],z(ind2,ind1),'linestyle','none');
    end
end
set(hax0,'xlim',[min(vx) max(vx)],'ylim',[min(vy) max(vy)],'xtick',x,'ytick',y);
if nargout>0
    varargout{1} = hax0;
    varargout{2} = hpat;
end

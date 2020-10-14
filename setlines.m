function setlines(varargin)

% set linecolor and linestyle automatically
% max number of lines for a graph is 28

linestyle{1}='-';
linestyle{2}='--';
linestyle{3}=':';
linestyle{4}='-.';

if nargin==0
    h=get(gca,'children');
    h=flipdim(map2colvec(h),1);
elseif nargin==1
    h=get(varargin{1},'children');
end
N=size(h,1);
% load COLORS_line;
COLORS_line = generate_colormap(N);
clrs=COLORS_line;
Ncl=size(clrs,1);

if N~=0
    for ind2=1:N
        set(h(ind2),'color',clrs(mod(ind2-1,Ncl)+1,:),'linestyle',linestyle{ceil(ind2/Ncl)})
    end
end
end

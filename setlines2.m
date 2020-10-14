function setlines(varargin)

% set colorstyle automatically
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
clrs=cell2mat(get(h,'color'));
N=0;
for ind1=2:size(clrs,1)
    if clrs(1,:)*clrs(ind1,:).'==clrs(1,:)*clrs(1,:).'
        N=ind1-1;
        break;
    end
end
if N~=0
    n=floor(size(clrs,1)/N);
    for ind2=1:n;
        set(h(1+(ind2-1)*N:ind2*N),'linestyle',linestyle{ind2})
    end
    ind2=ind2+1;
    set(h(1+(ind2-1)*N:end),'linestyle',linestyle{n+1});
end
end

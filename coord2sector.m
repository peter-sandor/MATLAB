function varargout = coord2sector(filename,center,Nsector,Nbin,norm_mode,varargin)

if nargin==6
    Ecoeff=varargin{1};
end
cutoff=0.003;
data=dlmread(filename);
N=size(data);
[a,b,c]=fileparts(filename);
if ~isempty(a)
    a=[a '\'];
end
frags=unique(data(:,1));
frags(frags==0)=[];
Nfrag=length(frags);
dist=sort(sqrt((data(:,2)-center(1)).^2+(data(:,3)-center(2)).^2),'descend');
Rmax=dist(ceil(N(1)*cutoff));
if nargin==5
    Rbins=0:Rmax/(Nbin-1):Rmax;
elseif nargin==6
    Rbins=0:Ecoeff*Rmax^2/(Nbin-1):Ecoeff*Rmax^2;
end
Sbins=(0:2*pi/Nsector:2*pi*(1-1/(Nsector)));
out=zeros(length(frags),Nsector,Nbin);
for ind1=1:N(1)
    alpha=coord2angle(data(ind1,2:3)-center);
    radius=sqrt((data(ind1,2)-center(1)).^2+(data(ind1,3)-center(2)).^2);
    sectind=sum(double(alpha>Sbins));
    if nargin==5
        radiusind=sum(double(radius>Rbins));
    elseif nargin==6
        radiusind=sum(double(Ecoeff*radius^2>Rbins));
    end
    if sectind==0
        sectind=1;
    end
    out(data(ind1,1),sectind,radiusind)=out(data(ind1,1),sectind,radiusind)+1;
end

if strcmp(norm_mode,'total')
	varargout{1}=out./permute(extend(extend(sum(sum(out,3),1),Nbin),Nfrag),[3 1 2]);
elseif strcmp(norm_mode,'frag')
	varargout{1}=out./extend(sum(out,3),Nbin);
elseif strcmp(norm_mode,'none')
	varargout{1}=out;
end

if nargin==5
    varargout{2}=Rmax;
elseif nargin==6
    varargout{2}=Ecoeff*Rmax^2;
end

end
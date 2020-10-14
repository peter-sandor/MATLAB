function autocolorP(varargin)

% set linecolor automatically
if nargin==0
    h=get(gca,'children');
    h=flipdim(map2colvec(h),1);
elseif nargin==1
    h=get(varargin{1},'children');
end
N=size(h,1);
clrs0=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1;]; % start from red, go through yellow, green, teal, blue and purple then back to red
clrs1=clrs0;
Ncol=size(clrs0,1);
Ndiv=floor(N/Ncol);
% clrs=Bezier_curve(clrs0,N);
if N<=Ncol
    clrs=clrs0(1:N,:);
else
    Nmod=mod(N,Ncol);
    clrs=zeros([N 3]);
    ind1=1;
    for ind2=1:Ncol
        if ind2>Nmod
            flagmod=0;
        else
            flagmod=1;
        end
        clrs(ind1,:)=clrs0(ind2,:);
        ind1=ind1+1;
        for ind3=1:Ndiv-1+flagmod
            clrs(ind1,:)=ind3*((clrs1(2,:)-clrs1(1,:)))/(Ndiv+flagmod)+clrs1(1,:);
            ind1=ind1+1;
        end
        clrs1=circshift(clrs0,-ind2,1);
	end
end

if N~=0
    for ind2=1:N;
        set(h(ind2),'color',clrs(ind2,:),'linestyle','-')
    end
end
end

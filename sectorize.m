function varargout = sectorize(imgin,center,Nphi,Nrad,varargin)

% This code sorts the values of a 2D rectangular array into radial and
% angular bins.

% input parameters:
% imgin: 2D rectangular array to be sorted
% center: coordinates [x y] a.k.a. [row column]
% Nphi: number of angular sectors
% Nrad: number of radial sectors
% interpfactor: each pixel will be subdivided into 2^interpfactor pixels.
%           0 means no subdivision
% flagR: if set to 1, radial lineouts will be multiplied by radius, for 0
%           they won't.
% Ecoeff: if set to other than 0, width of radial bins will be not equal,
%           but linearly increasing with distance from center. In other
%           words, bins will have equal size in a quadratic space. This option can
%           be used for transforming VMI distributions directly from momentum to
%           energy space, where the transformation is defined as: E=Ecoeff*pixel^2
%
% varargin{1:4} = interpfactor,flagR,Ecoeff,radmode

if nargin==4 % same behaviour as radial_lineoutVer3P
    interpfactor=0;
    flagR=1;
    Ecoeff=0;
    radmode=1; %'minRad' or 'maxRad';
else
    interpfactor=varargin{1};
    flagR=varargin{2};
    Ecoeff=varargin{3};
    if strcmp(varargin{4},'maxrad')==1
        radmode=1;
    else
        radmode=0;
    end
end
N0=size(imgin);
imip=interp2(imgin,interpfactor);
N=size(imip);
for ind1=1:interpfactor
    center=2*center-[1 1];
end

if radmode==0;
    Rmax=min([center(1) N(2)-center(1) center(2) N(1)-center(2)]);
elseif radmode==1
    corners=[1 1; 1 N(1); N(2) 1; N(2) N(1)];
    Rmax=max(sqrt(sum((permute(extend(center,4),[2 1])-corners).^2,2)));
end
% if Nrad==0
%     Nrad=round(Rmax);
% end
samplingfactor=Nrad/Rmax;
if Ecoeff==0
    Rbins=map2colvec(0:Rmax/(Nrad-1):Rmax)/samplingfactor; % calculate bins in momentum space
else
    Ecoeff=Ecoeff/samplingfactor;
    Rbins=map2colvec(0:Ecoeff*Rmax^2/(Nrad-1):Ecoeff*Rmax^2); % calculate bins in energy space
end
Sbins=map2colvec(0:2*pi/Nphi:2*pi*(1-1/(Nphi)));
out=zeros([Nrad Nphi]);
sect_regist=zeros([Nrad Nphi]);
for ind1=1:N(1)
    for ind2=1:N(2)
         % for a given pixel (ind1,ind2) calculate the angle with respect
         % to the edge of the image (alpha), and calculate the distance from the
         % center (radius)
        alpha = angle(+i*(ind2-center(1))+(ind1-center(2)))+pi;
        radius=sqrt((ind2-center(1)).^2+(ind1-center(2)).^2);
        if Ecoeff~=0
            radius=Ecoeff*radius^2;
        end
        if ~isempty(alpha)
            % given alpha and radius, figure out which bin should be added
            % to contribution to --> (radiusind,sectind)
            sectind=sum(double(alpha>Sbins));
            radiusind=sum(double(radius>Rbins));
            if radiusind==0
                radiusind=1;
            end
            if sectind==0
                sectind=1;
            end
        else
            sectind=1:length(Sbins);
            radiusind=1;
        end
        out(radiusind,sectind)=out(radiusind,sectind)+imip(ind1,ind2); % add the contribution (=pixel value) to the appropriate bin
        sect_regist(radiusind,sectind)=sect_regist(radiusind,sectind)+1; % also keep a record of how many pixels contributed to the given bin
    end
end
sect_regist = (sect_regist==0) + sect_regist; % make sure we don't divide by zero
if flagR==1
    for indS = 1 : Nphi
        out(:,indS) = out(:,indS) ./sect_regist(:,indS).*Rbins;
    end
else
    out=out./sect_regist; % normalize by the number of pixels contributing to each bin
end
varargout{1}=out;
varargout{2}=Rbins;
varargout{3}=Sbins;
end
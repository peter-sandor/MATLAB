function varargout = coord_binner(coords,center,Nphi,Nrad,Ecoeff)

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

flagR=1;
Nc=size(coords,1);
Rmax=max(sqrt((coords(:,1)-center(1)).^2+(coords(:,2)-center(2)).^2));
samplingfactor=Nrad/Rmax;
Ecoeff=Ecoeff/samplingfactor^2;
if Ecoeff==0
    Rbins=map2colvec(0:Rmax/(Nrad-1):Rmax)/samplingfactor; % calculate bins in momentum space
else
    Rbins=map2colvec(0:Ecoeff*Rmax^2/(Nrad-1):Ecoeff*Rmax^2); % calculate bins in energy space
end
Sbins=map2colvec(0:2*pi/Nphi:2*pi*(1-1/(Nphi)));
out=zeros([Nrad Nphi]);
sect_regist=zeros([Nrad Nphi]);
for ind1=1:Nc
    alpha = angle(+i*(coords(ind1,1)-center(2))+(coords(ind1,2)-center(1)))+pi;
    radius=sqrt((coords(ind1,1)-center(2)).^2+(coords(ind1,2)-center(1)).^2);
    if Ecoeff~=0
        radius=Ecoeff*radius^2;
    end
    if ~isempty(alpha)
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
%         disp([num2str(sectind) ' ' num2str(radiusind)])
    out(radiusind,sectind)=out(radiusind,sectind)+1;
end
if flagR==1
    for indS = 1 : Nphi
        out(:,indS) = out(:,indS).*Rbins;
    end
end
varargout{1}=out;
varargout{2}=Rbins;
varargout{3}=Sbins;
end
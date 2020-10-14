function varargout = peak_props(varargin)
% takes a single-peaked distribution which is background subtracted and normalized to max value then calculates the mean, the sigma, and the FWHM.
% first argument should be an N-by-2 matrix
% first column: x-values
% second column: y-values

thrs=0.03;
datain=varargin{1};
N=size(datain,1);
meanx=sum(datain(:,1).*datain(:,2))/sum(datain(:,2));
mean_ind=round(sum(map2colvec(1:N).*datain(:,2))/sum(datain(:,2)));
sigma=sqrt(sum((datain(:,1)-meanx).^2.*datain(:,2))/sum(datain(:,2)));
% sigmay=0.01;
ind_divide=round(mean(vec2ind(datain(:,2)==max(datain(:,2)))));
point1=round(mean([min(vec2ind(datain(1:ind_divide,2)>0.5)) max(vec2ind(datain(1:ind_divide,2)<0.5))]));
point2=ind_divide + round(mean([max(vec2ind(datain(ind_divide+1:end,2)>0.5)) min(vec2ind(datain(ind_divide+1:end,2)<0.5))]));
% point2=sum(datain(:,1).*exp(-(datain(:,2)-0.5).^2/sigmay^2).*(datain(:,1)>=meanx))/sum(exp(-(datain(:,2)-0.5).^2/sigmay^2).*(datain(:,1)>=meanx)); %abscissa points for the half-maximum values
% point1=sum(datain(:,1).*exp(-(datain(:,2)-0.5).^2/sigmay^2).*(datain(:,1)<=meanx))/sum(exp(-(datain(:,2)-0.5).^2/sigmay^2).*(datain(:,1)<=meanx));
indvec1=unique((datain(:,2)>thrs).*map2colvec(1:N));
indvec1(indvec1==0)=[];
point3=max(1,round(4*(min(indvec1)-mean_ind)+mean_ind));
point4=min(N,round(4*(max(indvec1)-mean_ind)+mean_ind));
FWHM=datain(point2,1)-datain(point1,1);

varargout{1}=meanx;
varargout{2}=sigma;
varargout{3}=FWHM;
varargout{4}=[point1 point2]; % abscissa points for FWHM values
varargout{5}=[point3 point4]; % abscissa points for region where values are nonzero
end
function varargout = InvAbelImg(imagein,center,varargin)

% This code takes an image, folds it to one quadrant, and calculates inverse-Abel transform columnwise using the
% 1D algorithm 'invAbel5'. By default, it is assumed that the laser is
% polarized in 'x' direction (along a row vector)
% Input arguments:
%   imagein: 2D array (VMI image)
%   center: coordinates for the image center: [x y] a.k.a. [row column] (can be noninteger)
%   optional argument #1: if the character 'y' is passed, it is assumed that the laser is polarized along 'y' axis,and the folded image
%       is transposed
% Output arguments:
%   output argument #2: Abel-inverted image, number of shots (image multiplied by rho)
%   output argument #1: Abel-inverted image, densities
%   output argument #3: radial lineouts for the different sectors

folded=foldimage(imagein,center);
if nargin==3 && strcmp(varargin{1},'y')
    folded=folded.';
    center=flipdim(center,2);
end
N=size(folded);
for ind3=1:N(2)
	[rho,f0(:,ind3)] = invAbel5(1:N(1),map2rowvec(squeeze(folded(:,ind3))),N(1));
end
% unfold image, calculate radial lineouts
f1=f0.*(map2colvec(0:(size(f0,1)-1))*ones([1 size(f0,2)]));
M=2*(N-1)+1;
% varargout{3}=sectorize(unfoldimage(f1,4),N,10);
varargout{1}=unfoldimage(f0,4);
varargout{2}=unfoldimage(f1,4);
end
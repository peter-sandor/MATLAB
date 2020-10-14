function varargout = AbelImg(imagein,center,varargin)

% This code takes an image, folds it to one quadrant, and calculates forward-Abel transform columnwise using the
% 1D algorithm 'Abel'. By default, it is assumed that the laser is
% polarized in 'x' direction (along a row vector)
% Input arguments:
%   imagein: 2D array (VMI image)
%   center: coordinates for the image center: [x y] a.k.a. [row column] (can be noninteger)
%   optional argument #1: if the character 'y' is passed, it is assumed that the laser is polarized along 'y' axis,and the folded image
%       is transposed
% Output arguments:
%   output argument #3: forward-projected image, number of shots (image multiplied by rho)
%   output argument #2: forward-projected image, densities
%   output argument #1: radial lineouts for the different sectors

folded=foldimage(imagein,center);
if nargin==3 && strcmp(varargin{1},'y')
    folded=folded.';
    center=flipdim(center,2);
end
N=size(folded);
for ind3=1:N(2)
	[rho,F0(:,ind3)] = Abel(1:N(1),map2rowvec(squeeze(folded(:,ind3))),N(1));
    disp(['line: ' num2str(ind3) '/' num2str(N(2))]);
end
% unfold image, calculate radial lineouts
% f1=f0.*(map2colvec(0:(size(f0,1)-1))*ones([1 size(f0,2)]));
M=2*(N-1)+1;
% varargout{1}=radial_lineoutVer3P(unfoldimage(f1,4),N,10);
varargout{1}=unfoldimage(F0,4);
% varargout{2}=unfoldimage(f1,4);
end
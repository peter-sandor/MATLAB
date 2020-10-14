function yout = smooth1D(yin,sigma)

% subroutine to Gaussian smooth a 1D array yin
% -----------------
% input parameters:
% yin       - array containing data to be smoothed
% sigma     - HWHM of gaussian used for the convolution (given in indices)
% -----------------
% output parameters:
% yout      - array of the smoothed data
yin=map2colvec(yin);
N=length(yin);
yout=zeros(size(yin));
vec=map2colvec(1:N);

for k=1:N
    func_conv = 1/(sqrt(2*pi)*sigma)*exp(-(vec-k).^2/(2*sigma^2));
    factor=sum(func_conv); % renormalizing the gaussian to remove the artifact at the ends of the smoothed data array
	yout(k) = sum(func_conv.*yin)/factor;
end
end
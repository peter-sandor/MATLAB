function varargout = clean_up_spectrum(data_in,varargin)
% Apply supergaussian filter to clean up spectrum
A = varargin{1};
if min(size(A)) > 1 % provide reference curve from which parameters for the supergaussian filter can be determined
%     x0 = peak_props([A(:,1) A(:,2).^2]);
    index_nonzero = logical(A(:,2)>0);
    x0_SG = (min(A(index_nonzero,1)) + max(A(index_nonzero,1)))/2;
    order_SG = 8;
    sigma_SG = fit_SG_sigma(A(:,1),index_nonzero,order_SG,x0_SG); % determine width of supergaussian filter
elseif min(size(A)) == 1 % provide parameters for supergaussian filter directly
    order_SG = A(1);
    sigma_SG = A(2);
    x0_SG = A(3);
end
window_SG = map2colvec(exp(-((data_in(:,1)-x0_SG).^2/sigma_SG^2).^order_SG));
data_out(:,1) = data_in(:,1);
data_out(:,2) = data_in(:,2).*window_SG;
varargout{1} = data_out;
varargout{2} = [order_SG sigma_SG x0_SG];
end
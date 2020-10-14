function FIT = linear_fit_to_data(x,y,varargin)

Nx = length(x);
fit_guess(2) = (y(end)-y(1))/(x(end)-x(1));
if abs(x(end)) > abs(x(1))
    fit_guess(1) = -fit_guess(2)*x(1)+y(1);
else
    fit_guess(1) = -fit_guess(2)*x(end)+y(end);
end
FIT = ezfit(x,y,['y = A + B*x; A = ' num2str(fit_guess(1)) '; B = ' num2str(fit_guess(2))]);
if nargin == 3
    stats.sigma_y = varargin{1};
    stats.F = [ones([Nx 1]) map2colvec(x)];
    F_T = stats.F';
    stats.W = diag(1./stats.sigma_y.^2);
    stats.R = F_T*stats.W*stats.F;
    % R = zeros(2);
    % R(1,2) = sum(y./map2rowvec(sigma_y).^2);
    % R(2,1) = R(1,2);
    % R(2,2) = sum(y.^2./map2rowvec(sigma_y).^2);
    % R(1,1) = sum(1./map2rowvec(sigma_y).^2);
    stats.R_inv = inv(stats.R);
    stats.chi2 = sum((map2colvec(FIT.m(1) + FIT.m(2)*x) - map2colvec(y)).^2);
    stats.B = stats.chi2/(Nx-length(FIT.m))*stats.R_inv;
    stats.coeffs = [FIT.m(1) sqrt(stats.B(1,1)); FIT.m(2) sqrt(stats.B(2,2))];
    [stats.y2, stats.sigma_y2] = hist_from_linear_fit(x,stats.coeffs);
    FIT.stats = stats;
end
end
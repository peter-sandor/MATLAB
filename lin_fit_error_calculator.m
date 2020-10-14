function FIT_error = lin_fit_error_calculator( x_data, y_data, axis_cross, slope, varargin )
%lin_fit_error_calculator
% determine the the error for "a" and "b" for a linear fit (a*x+b)
% assumption: x values are precise(no error)
%calculation is based on linear chi square fitting

a = axis_cross;
b = slope;

if nargin == 5
    sigma_i = varargin{1};
    
    s = sum(sigma_i.^2);
    s_x = sum(x_data./(sigma_i.^2));
    s_xx = sum((x_data.^2)./(sigma_i.^2));
    delta = s*s_xx-(s_x)^2;
    sigma_a_square = abs(s_xx/delta);
    sigma_b_square = abs(s/delta);
    
else
    sigma_i = ones(length(x_data),1);

    s = sum(sigma_i.^2);
    s_x = sum(x_data./(sigma_i.^2));
    s_xx = sum((x_data.^2)./(sigma_i.^2));
    delta = s*s_xx-(s_x)^2;
    chi_square = sum(((y_data-(a+x_data*b)).^2)./sigma_i);
    sigma_square = chi_square/(length(x_data)-2);
    sigma_a_square = (s_xx/delta)*sigma_square;
    sigma_b_square = (s/delta)*sigma_square;

end

% axis_cr_err = sqrt(sigma_a_square);
% slope_err = sqrt(sigma_b_square);

error.axis_cr =sqrt(sigma_a_square);
error.slope =sqrt(sigma_b_square);

FIT_error.error = error;

end


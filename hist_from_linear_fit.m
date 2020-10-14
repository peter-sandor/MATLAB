function [y,sigma_y] = hist_from_linear_fit(x,coeffs)
% This code calculates a histogram from the uncertainty of the values of 'y'
% for the linear fit y = a0 + a1*x. The input parameter contains the fitted
% parameters and their standard deviations: coeffs = [a0 sigma_a0; a1 sigma_a1]

N = 1000;
Nbins = 50;
Nx = length(x);
y = coeffs(1,1) + coeffs(2,1)*map2colvec(x);
% If X is a stochastic variable with standard normal distribution (E(X)=0,
% D2(X)=1), then Y=a*X+b is a stochastic variable with normal distribution but with
% expectation value 'b' (E(Y) = b) and variance 'a^2' (D2(Y) = a^2).
for ind1 = 1:N
    lines_rand(:,ind1) = coeffs(1,1)+randn*coeffs(1,2) + (coeffs(2,1)+randn*coeffs(2,2))*x;
end
sigma_y = std(lines_rand,[],2);
if 0
    ind_check = round(Nx/2);
    vec_bins = min(lines_rand(ind_check,:)):(max(lines_rand(ind_check,:))-min(lines_rand(ind_check,:)))/(Nbins-1):max(lines_rand(ind_check,:));
    counts = hist(lines_rand(ind_check,:),vec_bins);
    FIT1 = ezfit(vec_bins,counts,['A*exp(-(x-B)^2/(2*C^2));A=' num2str(max(counts)) ';B=' num2str(y(ind_check)) ';C=' num2str(sigma_y(ind_check))]);
    figure;bar(vec_bins,counts)
    showfit(FIT1)
end
end
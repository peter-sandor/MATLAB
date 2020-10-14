function varargout = expand_in_Legendre(x_in,y_in,maxorder)

prob_thrs = 0.05;
N = length(x_in);
x_center = (x_in(end) + x_in(1))/2;
x_span = x_center - x_in(1);
x = roundP((x_in-x_center)/x_span,6);
LegP = LegendreP(x,maxorder,'norm');
for ind0 = 1:maxorder+1
    y = y_in;
    y_expanded{ind0} = zeros(size(x));
    for ind1 = 1:ind0
        if ind1>1
            y = y - coeffs{ind0}(ind1-1)*squeeze(LegP(:,ind1-1,maxorder+1));
        end
        coeffs{ind0}(ind1) = trapz(x,y.*squeeze(LegP(:,ind1,maxorder+1)));
        y_expanded{ind0} = y_expanded{ind0} + coeffs{ind0}(ind1)*squeeze(LegP(:,ind1,maxorder+1));
    end
    chi2(ind0) = 1/(N-ind0)*sum((y_expanded{ind0} - y_in).^2);
    if ind0>2
        [prob1(ind0), Fstat] = F_test(N,ind0,ind0-1,chi2(ind0),chi2(ind0-1));
    else
        prob1(ind0) = 0;
    end
    if 0
        figure;plot(x_in,y_in,'k')
        hold on;plot(x_in,y_expanded{ind0},'r')
        title(num2str(['N_{basis} = ' num2str(ind0) '; \chi^2 = ' num2str(chi2(ind0)) ', ' num2str(prob1(ind0))]))
    end
end

ind_min_chi2 = vec2ind(chi2==min(chi2));
% ind_select = min(vec2ind(prob1>prob_thrs))-1;
ind_select = ind_min_chi2;
% disp(['probs = ' num2str(prob1) '; ind_select = ' num2str(ind_select) '; ind_min_chi2 = ' num2str(ind_min_chi2)]);
varargout{1} = coeffs{ind_select};
varargout{2} = y_expanded{ind_select};
varargout{3} = x;
end
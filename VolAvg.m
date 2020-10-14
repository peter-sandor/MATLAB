function varargout = VolAvg(intensity)
% 'I0' is the fractional intensity compared to the peak. Can not be 0.
% 'frac_vol' is the fractional volume that 'sees' a higher intensity than
% 'Int'.
% dV is the differential volume element that sees 'Int'.
% Source: Marc Hertlein's thesis, Appendix C.3
frac_vol=pi*(4/3*(1./intensity-1).^(1/2)+2/9*(1./intensity-1).^(3/2)-4/3*atan((1./intensity-1).^(1/2)));
frac_vol=frac_vol/sum(frac_vol);
dV=(2./(3*intensity)+1./(3*intensity.^2)).*sqrt(1./intensity-1);
% dV=dV/sum(dV);

varargout{1}=dV;
varargout{2}=frac_vol;
end
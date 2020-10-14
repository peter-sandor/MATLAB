function field_out = apply_spectral_mask(omega,field_in,M)
% omega_in:   angular frequencies
% field_in:   field in spectral domain (can be complex)
% M:        string describing the field mask in the spectral domain as a function of frequency. Use
%           'omega' to refer to angular frequency as the independent variable and use omega0 to
%           denote the central angular frequency (to be determined in this code).
% field_out: electric field vs omega with spectral mask applied

omega0 = sum(omega.*abs(field_in).^2)/sum(abs(field_in).^2);
eval(['mask=' M ';']);
field_out = field_in.*map2colvec(mask);
end
function intout = chirp2int(pulsedur,chirp)
% calculates relative intensity compared to transform limited case for
% different applied chirps to a Gaussian pulse. Source:
% [www.rp-photonics.com/chromatic_dispersion.html]
% Input parameters:
%   pulsedur    - transform limited intensity FWHM pulse duration in [fs]
%   chirp       - applied chirp in [fs^2] (--> d^2(Phi)/d(omega^2), Phi(omega)=...+1/2*d^2(Phi)/d(omega^2)*(omega-omega0)^2 )
% Output parameters:
%   intout      - relative intensity (dimensionless)

intout=1./sqrt(1+(4*log(2)*chirp./pulsedur.^2).^2);
end
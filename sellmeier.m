function varargout = sellmeier(lambda,material)

% This code calculates the refractive index of various materials as a
% function of wavelength.
% input parameters: 'lambda' -  wavelength in [micron], Nx1 or 1xN array
%                   'material'  -  material identifier
% output:           'n'  -  refractive indices
% available material identifiers:
% 'BK7' - borosilicate glass
% 'FS' - Fused Silica
% 'SF10'
% 'BBO_o' - Beta-Barium Borate, ordinary polarization
% 'BBO_e' - Beta-Barium Borate, extraordinary polarization
% 'TeO2' - Tellurium Dioxide
% 'GaP' - Gallium Phosphate
% Calcite_o - crystalline calcite ordinary polarization
% Calcite_e - crystalline calcite extraordinary polarization
% 'PbMoO4'
% 'Quartz_o'
% 'Quartz_e
% 'Sapphire'

if strcmp(material,'FS')
    lrange=[0.21 3.71];
    n0=1;
    coeffs=[0.6961663 0.0684043^2; 0.4079426 0.1162414^2; 0.8974794 9.896161^2];
% source: http://www.rp-photonics.com/sellmeier_formula.html
% works between 0.21 and 3.71 um
elseif strcmp(material,'BK7');
    lrange=[0.249 2.325];
    n0=1;
    coeffs=[1.03961212 0.00600069867; 0.231792344 0.0200179144; 1.01046945 103.560653];
% source: refractiveindex.info
% works between 0.4 to 1 um
elseif strcmp(material,'BBO_e');
    lrange=[0.40 1.10];
    n0=1;
    coeffs=[0.5525 0; 0.8205 0.1249^2; 0.423 10^2];
elseif strcmp(material,'BBO_o');
    lrange=[0.40 1.10];
    n0=1;
    coeffs=[0.7126 0; 1.0279 0.1338^2; 1.535 10^2];
elseif strcmp(material,'TeO2');
    lrange=[0.40 1.00];
    n0=1;
    coeffs=[2.584 0.1342^2; 1.157 0.2638^2];
% source: refractiveindex.info
% works between 0.4 to 1 um
elseif strcmp(material,'PbMoO4');
    lrange=[0.44 1.08];
    n0=1;
    coeffs=[3.54642 0.18518^2; 0.58270 0.033764^2];
% source: refractiveindex.info
% works between 0.44 to 1.08 um
elseif strcmp(material,'SF10');
    lrange=[0.38 2.5];
    n0=1;
    coeffs=[1.62153902 0.0122241457; 0.256287842 0.0595736775; 1.64447552 147.468793];
% source: refractiveindex.info
% works between 0.38 to 2.5 um
elseif strcmp(material,'GaP');
    lrange=[0.5 0.825];
    n0=1;
    coeffs=[7.2723 0.06233; 0.86832 0.15127];
% source: refractiveindex.info, tabular data. Coefficients obtained from manual fit.
% works between 0.5 to 0.825 um
elseif strcmp(material,'Calcite_o');
    lrange=[0.204 2.172];
    n0=1;
    coeffs=[0.73358749 0; 0.96464345 1.94325203e-2; 1.82831454 120];
% source: refractiveindex.info, Ghosh 1999.
% works between 0.204 to 2.172 um
elseif strcmp(material,'Calcite_e');
    lrange=[0.204 2.172];
    n0=1;
    coeffs=[0.35859695 0; 0.82427830 1.06689543e-2; 0.14429128 120];
% source: refractiveindex.info, Ghosh 1999.
% works between 0.204 to 2.172 um
elseif strcmp(material,'Quartz_o');
    lrange=[0.204 2.1];
    n0=1.28604141;
    coeffs=[1.07044083 1.00585997*1e-2; 1.10202242 100;];
% source: Gorachand Gosh: Dispersion-equation coefficients for the refractive index and
% birefringence of calcite and quartz crystals (Optics Communications 163 1999 95–102)
% works between 0.204 to 2.1 um
elseif strcmp(material,'Quartz_e');
    lrange=[0.204 2.1];
    n0=1.28851804;
    coeffs=[1.09509924 1.02101864*1e-2; 1.15662475 100;];
% source: Gorachand Gosh: Dispersion-equation coefficients for the refractive index and
% birefringence of calcite and quartz crystals (Optics Communications 163 1999 95–102)
% works between 0.204 to 2.1 um
elseif strcmp(material,'Sapphire_o');
    lrange=[0.2 5];
    n0=1;
    coeffs=[1.4313493 0.0726631^2; 0.65054713 0.1193242^2; 5.3414021 18.028251^2];
% source: https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire
elseif strcmp(material,'Air');
    lrange=[0.23 1.69];
    n0=1;
    coeffs=[0.05792105/238.0185 1/238.0185; 0.00167917/57.362 1/57.362];
% Standard air: dry air at 15 °C, 101.325 kPa and with 450 ppm CO2 content.
% P. E. Ciddor. Refractive index of air: new equations for the visible and near infrared, Appl. Optics 35, 1566-1573 (1996)
else
    disp('Material is not in the catalog.');
    return;
end

if all((lambda>=min(lrange)).*(lambda<=max(lrange)))
    n=n0;
    for ind1=1:size(coeffs,1)
        n = n + coeffs(ind1,1) * lambda.^2 ./(lambda.^2 - coeffs(ind1,2));
    end
    n=sqrt(n);
    varargout{2}=lambda;
    varargout{1}=n;
else
    varagout{1}=[];
    varargout{2}=[];
    disp(['Wavelength values out of range! (' num2str(lrange(1)) ' - ' num2str(lrange(2)) ' micron)']);
end
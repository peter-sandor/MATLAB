function evaluate_FE_results(varargin)
if nargin == 1
    folder_name = varargin{1};
else
    folder_name = [];
end
temp = load([folder_name 'data.mat']);
temp_fieldname = fieldnames(temp);
data = getfield(temp,temp_fieldname{1});

constants_si;
peak_intensity = data.result.PeakIntensity*1e4; % [W/m^2]
E_cutoff = data.result.E_cutoff*C_SI.qE; % [J]
lambda = data.int_calib.lambda_nm*1e-9; % [m]
WorkFunction = data.int_calib.WorkFunction_eV; % [eV]
FE_2 = 1*mean(data.result.field_enh);
field = sqrt(2*peak_intensity/C_SI.eps0/C_SI.c); % [V/m]
omega = 2*pi*C_SI.c/lambda; % [Hz]

% Evaluate field_enhancement simply from the slope of the E_cutoff vs
% Intensity curve.
FIT1 = linear_fit_to_data(peak_intensity,E_cutoff);
E_cut_slope = FIT1.m(2);
FE = sqrt(E_cut_slope*2*C_SI.mE*C_SI.eps0*C_SI.c*omega^2/C_SI.qE^2/10.007);

Up = C_SI.qE^2/(4*C_SI.mE)*(FE.*field).^2/omega^2/C_SI.qE; % [eV]
E_cutoff_calc = 10.007*Up + 0.538*WorkFunction; % [eV]
keldysh = omega*sqrt(2*C_SI.mE*WorkFunction*C_SI.qE)./(C_SI.qE*FE.*field);

disp([num2str(FE) ' ' num2str(FE_2)]);

figure;
subplot(211);hold on;
plot(peak_intensity*1e-4,data.result.keldysh(data.result.keldysh~=0),'ko');
plot(peak_intensity*1e-4,keldysh,'r+')
xlabel('Peak intensity [W/cm^2]');
ylabel('\gamma_{Keldysh}');
title(['FE = ' num2str(FE)]);
legend('read from data save to disk','recalculated from FE from slope')

subplot(212);hold on;
plot(peak_intensity*1e-4,data.result.E_cutoff,'ko');
plot(peak_intensity*1e-4,E_cutoff_calc,'r+')
xlabel('Peak intensity [W/cm^2]')
ylabel('E_{cutoff} [eV]')
% showfit(FIT1)
end
function [energy,spectrum] = Kaesdorf_TOF_delay2energy(data,t0,calibfile)

data_calib = load(calibfile);
energy_calib = data_calib(:,1); % [eV]
delay_calib = data_calib(:,2)*1000; % [ns]
dtdE = num_diff(energy_calib,delay_calib);

delay0 = data(:,1) + t0; % [ns]
counts0 = data(:,2);
energy = map2colvec(interp1(delay_calib,energy_calib,delay0,'linear','extrap'));
dtdE_ip = interp1(energy_calib,dtdE,energy,'linear','extrap');
spectrum = map2colvec(-dtdE_ip.*counts0);
[energy,ind_sort] = sort(energy);
spectrum = spectrum(ind_sort);
end

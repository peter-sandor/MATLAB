%% Create directory and read raw data
purge;

if exist('params_for_analyze_FE_itx.m')==2
    params_for_analyze_FE_itx;
end

if 1
    E_fit3 = E_fit2;
end

if exist('params_for_analyze_FE_itx.m')==0
    str_files0 = 'spectrum*.itx';
    str_files1 = 'spectrum_*.itx';
    folder_to_save = 'figures4\';

    % Set 'flag_plateau' to one if you want to calculate 'E_high' (and hence the cutoff energy: E_cutoff=E_high-E0)
    % from the intersect between the fits to the plateau and high-energy tail of the spectrum.
    % Otherwise, 'E_high' is calculated from the energy where the
    % spectrum drops to '1/peak_part' value (see above).
    flag_plateau = 1;
    % Set 'flag_E0_from_fit' to 1 if you want to calculate 'E0' from a fit to
    % the low-energy slope of the spectrum. Set it to 2 if you would like a
    % manually determined constant value (e.g. E0 = 3.0 eV).
    % Otherwise it will be calculated as E0 = E_maxyield - 2.3 eV,
    % where 'E_maxyield' is the energy where the yield is maximum.
    flag_E0_from_fit = 2;
    flag_save = 1; % Set this to 1 if you want to save the figures automatically.
end

if exist(folder_to_save)==0
    mkdir(folder_to_save);
end

files0 = dir(str_files0);
if 0
    if ~isempty(strfind(files0(1).name,' '))
        for ind1 = 1:length(files0)
            filename_out = space2uscore(files0(ind1).name,'copy');
            pad_name_w_zeros(filename_out,2,'rename');
        end
    end
    files1 = dir(str_files1);
else
    files1 = files0;
end
nfile = length(files1);  % number of files
fontsize = 16;
linewidth = 1.5;
constants_si;
% Read data
for ind1 = 1:nfile
    [data_raw{ind1} nrg{ind1} alpha{ind1}] = read_itx_multiaxes(files1(ind1).name);
end
%% Preprocess data
 % Write power values in the same order as they appear in .sle file!
 % Write zero for background trace.
 
if exist('params_for_analyze_FE_itx.m')==0
    power = [40:-2.5:25 21.74 20 17.5];
    ind_select0 = [1:10];
    focus_x_um = 3.1; %[um] Radius; based on Top Hat criterion (see Siegman: Lasers, chapter 17.1, section title: 'Aperture Transmission')
    focus_y_um = 3.1; %[um]
    lambda = 800*10^-9; % [m]
    RepetitionRate = 80*10^6; %[Hz]
    PulseLength_fs = 7; %[fs]
    peak_part = 1e4;   %  for defining the yield level which is considered as the high-energy tail of the spectrum.
end

omega = 2*pi*C_SI.c/lambda;
eVtoJoule=1.602*10^-19;
WorkFunction_eV = 5.1;
cikkfactor = 0.538;

if ~isempty(vec2ind(power == 0))
    ind_bkg = vec2ind(power==0);
    N_trace = size(read_itx_multiaxes(files1(ind_bkg).name));
    ind_tr_bkg = [1 N_trace(2)];
else
    ind_bkg = vec2ind(power==min(power));
    N_trace = size(read_itx_multiaxes(files1(ind_bkg).name));
    ind_tr_bkg = [round(N_trace(2)/2) N_trace(2)];
end

if ~isempty(ind_bkg)
    data_bg = [map2colvec(sum(data_raw{ind_bkg},1))];
else
    [counts_bkg,nrg_bkg,alpha_bkg] = read_itx_multiaxes('?.itx');
    data_bg = [map2colvec(sum(counts_bkg,1))];
end
bkg_level = mean(data_bg(ind_tr_bkg(1):ind_tr_bkg(2)));
bkg_std = std(data_bg(ind_tr_bkg(1):ind_tr_bkg(2)));

[power_sorted,ind_sort] = sort(power(ind_select0));
ind_select1 = ind_select0(ind_sort);
Nspec = length(ind_select1);

Nbin = 50;
binwidth = (max(data_bg)-min(data_bg))/(Nbin-1);
bin_edges0 = (min(data_bg)-binwidth/2:binwidth:max(data_bg)+binwidth/2);
bin_centers = (min(data_bg):binwidth:max(data_bg));
[counts,bin_edges] = histcounts(data_bg,bin_edges0);
figure;
subplot(211);hold on;
plot(data_bg,'ko');
line([ind_tr_bkg(1) ind_tr_bkg(1)],[min(data_bg) max(data_bg)],'color','r');
line([ind_tr_bkg(2) ind_tr_bkg(2)],[min(data_bg) max(data_bg)],'color','r');
title('background scan')
subplot(212);bar(bin_centers,counts)
xlabel('signal level')
ylabel('counts')
title('Histogram of background trace')

E_max = 0;
if exist('params_for_analyze_FE_itx.m')==0
    E_cut = 9.5; % [eV]; energy below which the algorithm looks for the sharp low-energy peak
end
for ind2 = 1:Nspec
    Edata_plot{ind2}(:,1) = nrg{ind_select1(ind2)};
    Edata_plot{ind2}(:,2) = sum(data_raw{ind_select1(ind2)},1);
    Edata1{ind2}(:,1) = Edata_plot{ind2}(:,1);
    Edata1{ind2}(:,2) = Edata_plot{ind2}(:,2) - bkg_level;
    totalCounts(ind2) = sum(Edata1{ind2}(:,2));
    if 1
        Edata_norm{ind2} = Edata1{ind2}(:,2)/max(Edata1{ind2}(:,2));
    else
        Edata_norm{ind2} = Edata1{ind2}(:,2)/totalCounts(ind2);
    end
    Edata_n_sigma{ind2} = sqrt(Edata_plot{ind2}(:,2) + bkg_std^2)/totalCounts(ind2);
    Edata_plot{ind2}(Edata_plot{ind2}(:,2)<=0,2) = NaN;
    Ekin = Edata1{ind2}(:,1);
    ind_maxyield(ind2) = vec2ind(Edata_norm{ind2}(1:value2index(Ekin,E_cut)) == max(Edata_norm{ind2}(1:value2index(Ekin,E_cut))));
    E_maxyield0(ind2) = mean(ind_maxyield(ind2));
    if max(Edata1{ind2}(:,1)) > E_max
        E_max = max(Edata1{ind2}(:,1));
    end
end
% E_maxyield0(12) = 10.3;

PEAKPOS = figure;
plot(totalCounts,E_maxyield0,'ko')
xlabel('Total Counts [cps]')
ylabel('Position of peak maxima [eV]')
setfigP;
if flag_save
    saveas(PEAKPOS,[folder_to_save 'E_maxyielditions.fig']);
end

SPEKTRUM = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for ind2 = 1:Nspec
    plot(Edata_plot{ind2}(:,1),Edata_plot{ind2}(:,2),'Linewidth',linewidth)
%    plot(Edata_plot{ind2}(:,1),Edata_norm{ind2},'Linewidth',linewidth)
end
set(gca,'Yscale','log')
% xlim([0 E_max]);
xlabel('E_{kin} [eV]')
ylabel('Counts/s [au]')
title('Spectrum over laser power')
legend(cellstr(num2str(map2colvec(power_sorted)))+ " mW")
setfigP;
if flag_save
    saveas(SPEKTRUM,[folder_to_save 'spectra.fig']);
end
%% calculate peak intensity
if exist('params_for_analyze_FE_itx.m')==0
    ind_fit_int = [3 7];
end

FIT_yield0 = ezfit(log(power_sorted(ind_fit_int(1):ind_fit_int(2))), log(totalCounts(ind_fit_int(1):ind_fit_int(2))),'y = A + B*x; A = 0; B = 1e0;');
photon_order = FIT_yield0.m(2);
factor1 = exp(FIT_yield0.m(1));

AverageEnergy_uW = power_sorted*10^3; %[uW]
AverageEnergy_W = AverageEnergy_uW * 10^-6; %[W]
PulseLength_s = PulseLength_fs * 10^-15; %[s]

% focus_x_um = 1*sqrt(photon_order); %[um] RADIUS!!!
% focus_y_um = 1*sqrt(photon_order); %[um]
focus_x_cm = focus_x_um * 10^-4; %[cm]
focus_y_cm = focus_y_um * 10^-4; %[cm]

FocusSize = pi*focus_x_cm*focus_y_cm;   % [cm^2]
PulseEnergy = AverageEnergy_W / RepetitionRate; % [J]
PeakPower = PulseEnergy / PulseLength_s; % [W]
PeakIntensity = PeakPower / FocusSize; % [W/cm^2]
 
ElectricField = sqrt(2 / (C_SI.c*C_SI.eps0) * PeakIntensity*10^4 ); %[V/m]
ElectricField_nm = ElectricField*10^-9; %[V/nm]
FIT_yield1 = ezfit(log(PeakIntensity(ind_fit_int)), log(totalCounts(ind_fit_int)),'y = A + B*x; A = -1e2; B = 1e0;');
YIELD = figure;
% plot(log(PeakIntensity(ind_fit2(1):ind_fit2(2))),log(y_expfit(ind_fit2(1):ind_fit2(2))),'k+')
loglog(PeakIntensity,totalCounts(1:Nspec),'k+')
hold on;
plot(PeakIntensity(ind_fit_int),exp(FIT_yield1.m(1) + FIT_yield1.m(2)*log(PeakIntensity(ind_fit_int))),'r')
legend('Counts',sprintf('fit: y = %.3e *PeakInt^{%.4f}', exp(FIT_yield1.m(1)), FIT_yield1.m(2)));
xlabel('Peak Intensity [W/cm^2]');
ylabel('Electron Counts (a.u.)');
axis('tight')
set(gca,'Linewidth',linewidth,'FontSize',fontsize)
setfigP;
if flag_save
    saveas(YIELD,[folder_to_save 'yield_vs_intensity.fig']);
end
%% Main Loop with fits and Field Enhancement (FE) calculation

if 0
    temp = ginput(2*Nspec);
    E_fit_test = reshape(temp(:,1),[2 size(temp,1)/2]).'
end

if exist('params_for_analyze_FE_itx.m')==0
    E_fit0 = [
        8.5104    8.6714;...
        8.5409    8.7062;...
        8.5192    8.7062;...
        8.5148    8.7323;...
        8.5061    8.7149;...
        8.4887    8.7106;...
        8.2973    8.5540;...
        8.1058    8.3321;...
        7.4793    7.7795;...
        6.3046    6.5352]; % range for rising edge of narrow peak at low energy, [eV]

    E_fit1 = [
       10.2259   11.6766;...
        9.6854   10.4534;...
        9.7139   10.7379;...
        9.6285   11.4490;...
        9.6854   11.8188;...
        9.6285   11.5628;...
        8.9459   10.9654;...
        9.2587   11.3352;...
        9.3725   11.9894;...
        8.3770   10.1690]; % range for plateau, in [eV]

    E_fit2 = [
       11.5343   13.6392;...
       11.3637   13.2979;...
       12.1886   16.3699;...
       15.3744   18.3611;...
       15.5450   19.0153;...
       16.0002   20.3807;...
       18.0766   21.5753;...
       18.6455   22.8269;...
       19.2429   23.8225;...
       17.1664   23.1114]; % range for tail of spectrum, for intersection with plateau, in [eV]

    E_fit3 = [
       13.0988   17.2233;...
       13.5824   19.7833;...
       18.6740   24.1354;...
       21.3762   29.0279;...
       23.0545   30.5070;...
       24.6189   34.0342;...
       29.1985   36.8217;...
       31.6732   40.4342;...
       30.6208   39.3818;...
       29.3123   39.1827]; % range for tail of spectrum, for intersection with 1/peak_part, in [eV]

    ind_select2 = [3:Nspec];
end
SPECTRA_DETAILS=figure('units','normalized','outerposition',[0 0 1 1]);hax1=axes;  %to plot the baseline on this figure too, at the backround calc
hold on;

for ind1 = ind_select2
    Ekin = Edata1{ind1}(:,1);
    E_step = Ekin(2) - Ekin(1);
    hndl1(ind1) = plot(hax1,Ekin,Edata_norm{ind1},'marker','o');
    set(gca,'Yscale','log')
    ind_fit0(ind1,:) = value2index(Ekin,E_fit0(ind1,:));
    ind_fit1(ind1,:) = value2index(Ekin,E_fit1(ind1,:));
    ind_fit2(ind1,:) = value2index(Ekin,E_fit2(ind1,:));
    ind_fit3(ind1,:) = value2index(Ekin,E_fit3(ind1,:));
    ind_check1 = vec2ind(Edata_norm{ind1}(ind_fit3(ind1,1):ind_fit3(ind1,2))<=0);
    if ~isempty(ind_check1)
        ind_fit3(ind1,2) = ind_fit3(ind1,1) + min(ind_check1) - 2;
        E_fit3(ind1,2) = Ekin(ind_fit3(ind1,2));
    end
    [M I] = max(Edata_norm{ind1});
    norm_minmax(ind1,:) = [min(Edata_norm{ind1}) max(Edata_norm{ind1})];
    
    yield_cutlevel(ind1) = max(Edata_norm{ind1})/peak_part;
    log_cutlevel(ind1) = log(yield_cutlevel(ind1));
    std_cutlevel(ind1) = sqrt(yield_cutlevel(ind1));

    ind_zeros = vec2ind(Edata1{ind1}(ind_fit3(ind1,1):end,2)<=0);
    if ~isempty(ind_zeros)
        ind_zero(ind1) = min(ind_zeros) + ind_fit3(ind1,1);
    else
        ind_zero(ind1) = length(Ekin);
    end
    if ind_zero(ind1) < ind_fit3(ind1,2)
        ind_fit3(ind1,2) = ind_zero(ind1) - 2;
    end

	% Fit rising edge of electron spectrum with straight line and calculate
    % intersection with background
    x_linfit0 = Edata1{ind1}(ind_fit0(ind1,1):ind_fit0(ind1,2),1);
    y_linfit0 = log(Edata_norm{ind1}(ind_fit0(ind1,1):ind_fit0(ind1,2)));
    FIT_sp(ind1,1) = linear_fit_to_data(x_linfit0,y_linfit0,Edata_n_sigma{ind1}(ind_fit0(ind1,1):ind_fit0(ind1,2))./y_linfit0);
    intersect_linfit0(ind1) = (log(yield_cutlevel(ind1)) - FIT_sp(ind1,1).m(1))/FIT_sp(ind1,1).m(2);
    fit_errors(ind1,1) = lin_fit_error_calculator( x_linfit0, y_linfit0, FIT_sp(ind1,1).m(1), FIT_sp(ind1,1).m(2)); %calc error from results with chi square methode

    if flag_E0_from_fit == 1
        E0(ind1) = intersect_linfit0(ind1);
        E0_std(ind1) = sqrt(1/FIT_sp(ind1,1).m(2)^2*(FIT_sp(ind1,1).stats.coeffs(1,2)^2 + FIT_sp(ind1,1).m(1)^2/FIT_sp(ind1,1).m(2)^2*FIT_sp(ind1,1).stats.coeffs(2,2)^2));
    elseif flag_E0_from_fit == 2
        E0(ind1) = 3;
        E0_std(ind1) = 0;
    else
        E0(ind1) = E_maxyield0(ind1)-2.3;
        E0_std(ind1) = 0;
    end
    
    % Fit tail of the electron spectrum on log scale with a straight line
    % and calculate intersection with 1/peak_part signal level
    x_linfit3 = Edata1{ind1}(ind_fit3(ind1,1):ind_fit3(ind1,2),1);
    y_linfit3 = log(Edata_norm{ind1}(ind_fit3(ind1,1):ind_fit3(ind1,2)));
%     y_linfit2 = log(Edata1{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2),2));
    FIT_sp(ind1,2) = linear_fit_to_data(x_linfit3,y_linfit3,Edata_n_sigma{ind1}(ind_fit3(ind1,1):ind_fit3(ind1,2))./y_linfit3);
    intersect_linfit3(ind1) = (log(yield_cutlevel(ind1)) - FIT_sp(ind1,2).m(1))/FIT_sp(ind1,2).m(2);
    fit_errors(ind1,2) = lin_fit_error_calculator( x_linfit3, y_linfit3, FIT_sp(ind1,2).m(1), FIT_sp(ind1,2).m(2)); %calc error from results with chi square methode

    % Fit tail of the electron spectrum with exponential decay
    x_expfit = Edata1{ind1}(ind_fit3(ind1,1):length(Ekin),1);
    y_expfit = Edata_norm{ind1}(ind_fit3(ind1,1):length(Ekin));
    FIT_tail(ind1) = ezfit(x_expfit,y_expfit,['a*exp((x-' num2str(min(x_expfit)) ')*b)'], [max(y_expfit),-1]);
    intersect_expfit(ind1) = log(yield_cutlevel(ind1)/FIT_tail(ind1).m(1))/FIT_tail(ind1).m(2) + min(x_expfit);

    E_high(ind1) = (intersect_linfit3(ind1) + intersect_expfit(ind1))/2;
    E_high_std(ind1) = sqrt(1/FIT_sp(ind1,2).m(2)^2*(FIT_sp(ind1,2).stats.coeffs(1,2)^2 + FIT_sp(ind1,2).m(1)^2/FIT_sp(ind1,2).m(2)^2*FIT_sp(ind1,2).stats.coeffs(2,2)^2));
            
    switch flag_plateau
        case 0
            E_plateau(ind1) = 0;
            E_plateau_std(ind1) = 0;
            E_cutoff(ind1) = E_high(ind1) - E0(ind1);
            E_cutoff_std(ind1) = sqrt(E_high_std(ind1)^2 + E0_std(ind1)^2);
            E_cutoff_std_chi_sq(ind1) = 0; % for this to work, error calculation for the exponential fit needs to be implemented
        case 1
            % Fit plateau of the electron spectrum on log scale with a straight line
            % and calculate intersection with tail end
            x_linfit1 = Edata1{ind1}(ind_fit1(ind1,1):ind_fit1(ind1,2),1);
            y_linfit1 = log(Edata_norm{ind1}(ind_fit1(ind1,1):ind_fit1(ind1,2)));
            FIT_sp(ind1,3) = linear_fit_to_data(x_linfit1,y_linfit1,Edata_n_sigma{ind1}(ind_fit1(ind1,1):ind_fit1(ind1,2))./y_linfit1);
            fit_errors(ind1,3) = lin_fit_error_calculator( x_linfit1, y_linfit1, FIT_sp(ind1,3).m(1), FIT_sp(ind1,3).m(2)); %calc error from results with chi square methode

            x_linfit2 = Edata1{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2),1);
            y_linfit2 = log(Edata_norm{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2)));
            FIT_sp(ind1,4) = linear_fit_to_data(x_linfit2,y_linfit2,Edata_n_sigma{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2))./y_linfit2);
            fit_errors(ind1,4) = lin_fit_error_calculator( x_linfit2, y_linfit2, FIT_sp(ind1,4).m(1), FIT_sp(ind1,4).m(2)); %calc error from results with chi square methode

            E_plateau(ind1) = (FIT_sp(ind1,3).m(1) - FIT_sp(ind1,4).m(1))/(FIT_sp(ind1,4).m(2) - FIT_sp(ind1,3).m(2));
            E_plateau_std(ind1) = sqrt(1/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(FIT_sp(ind1,3).stats.coeffs(1,2)^2 + FIT_sp(ind1,4).stats.coeffs(1,2)^2 + (FIT_sp(ind1,3).m(1)-FIT_sp(ind1,4).m(1))^2/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(FIT_sp(ind1,3).stats.coeffs(2,2)^2 + FIT_sp(ind1,4).stats.coeffs(2,2)^2)));
            E_plateau_std_chi_sq(ind1) = sqrt(1/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(fit_errors(ind1,3).error.axis_cr^2 + fit_errors(ind1,4).error.axis_cr^2 + (FIT_sp(ind1,3).m(1)-FIT_sp(ind1,4).m(1))^2/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(fit_errors(ind1,3).error.slope^2 + fit_errors(ind1,4).error.slope^2)));
            E_cutoff(ind1) = E_plateau(ind1) - E0(ind1);
            E_cutoff_std(ind1) = sqrt(E_plateau_std(ind1)^2 + E0_std(ind1)^2);
            E_cutoff_std_chi_sq(ind1) = sqrt(E_plateau_std_chi_sq(ind1)^2 + E0_std(ind1)^2);
        case 2
            x_linfit2 = Edata1{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2),1);
            y_linfit2 = log(Edata_norm{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2)));
            FIT_sp(ind1,4) = linear_fit_to_data(x_linfit2,y_linfit2,Edata_n_sigma{ind1}(ind_fit2(ind1,1):ind_fit2(ind1,2))./y_linfit2);
            fit_errors(ind1,4) = lin_fit_error_calculator( x_linfit2, y_linfit2, FIT_sp(ind1,4).m(1), FIT_sp(ind1,4).m(2)); %calc error from results with chi square methode
            
            x_tail = Edata1{ind1}(ind_maxyield(ind1):ind_fit3(ind1,2),1);
            y_tail = log(Edata_norm{ind1}(ind_maxyield(ind1):ind_fit3(ind1,2)));
            [coeffs,y_Leg] = expand_in_Legendre(x_tail,y_tail,5);
            y_diff=num_diff(x_tail,y_Leg);
            ind_local_extr = value2index(y_diff,0);
            if 0
                figure;plot(x_tail,y_tail,'k')
                hold on;plot(x_tail,y_Leg,'r')
                figure;plot(x_tail,y_diff,'k')
            end
            if isempty(ind_local_extr)
                ind_local_extr = min(vec2ind(y_diff == max(y_diff)));
            end
            siglevel_local_extr(ind1) = y_tail(ind_local_extr);
            E_local_extr(ind1) = x_tail(ind_local_extr);
            E_plateau(ind1) = (siglevel_local_extr(ind1) - FIT_sp(ind1,4).m(1))/(FIT_sp(ind1,4).m(2));
            E_cutoff(ind1) = E_plateau(ind1) - E0(ind1);
            E_plateau_std(ind1) = 0; %sqrt(1/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(FIT_sp(ind1,3).stats.coeffs(1,2)^2 + FIT_sp(ind1,4).stats.coeffs(1,2)^2 + (FIT_sp(ind1,3).m(1)-FIT_sp(ind1,4).m(1))^2/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(FIT_sp(ind1,3).stats.coeffs(2,2)^2 + FIT_sp(ind1,4).stats.coeffs(2,2)^2)));
            E_plateau_std_chi_sq(ind1) = 0; %sqrt(1/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(fit_errors(ind1,3).error.axis_cr^2 + fit_errors(ind1,4).error.axis_cr^2 + (FIT_sp(ind1,3).m(1)-FIT_sp(ind1,4).m(1))^2/(FIT_sp(ind1,4).m(2)-FIT_sp(ind1,3).m(2))^2*(fit_errors(ind1,3).error.slope^2 + fit_errors(ind1,4).error.slope^2)));
            E_cutoff(ind1) = E_plateau(ind1) - E0(ind1);
            E_cutoff_std(ind1) =  E0_std(ind1); %sqrt(E_plateau_std(ind1)^2 + E0_std(ind1)^2);
            E_cutoff_std_chi_sq(ind1) =  E0_std(ind1); % sqrt(E_plateau_std_chi_sq(ind1)^2 + E0_std(ind1)^2);
    end
    
    x_plotax0 = E0(ind1):E_step/5:x_linfit0(end);
    plot(hax1,x_plotax0,exp(FIT_sp(ind1,1).m(2)*(E0(ind1):E_step/5:x_linfit0(end)) + FIT_sp(ind1,1).m(1)),'color','k')
    plot(hax1,x_plotax0,yield_cutlevel(ind1)*ones(size(x_plotax0)),'color','k')
    plot(hax1,E0(ind1),yield_cutlevel(ind1),'color','r','marker','+')
    
    plot(hax1,x_expfit,exp(FIT_sp(ind1,2).m(2)*x_expfit + FIT_sp(ind1,2).m(1)),'color','k')
    plot(hax1,x_expfit,yield_cutlevel(ind1)*ones(size(x_expfit)),'color','k')
    plot(hax1,intersect_linfit3(ind1),yield_cutlevel(ind1),'color','r','marker','+');
    if flag_plateau==1
        plot(hax1,x_linfit1,exp(FIT_sp(ind1,3).m(2)*x_linfit1 + FIT_sp(ind1,3).m(1)),'color','k')
        plot(hax1,x_linfit2,exp(FIT_sp(ind1,4).m(2)*x_linfit2 + FIT_sp(ind1,4).m(1)),'color','k')
        plot(hax1,E_plateau(ind1),exp(FIT_sp(ind1,3).m(2)*E_plateau(ind1) + FIT_sp(ind1,3).m(1)),'color','g','marker','+');
    elseif flag_plateau == 2
        plot(hax1,x_tail(ind_local_extr):(E_plateau(ind1) - x_tail(ind_local_extr)):E_plateau(ind1),exp(siglevel_local_extr(ind1))*ones([1 2]),'color','k','marker','none');
        plot(hax1,E_plateau(ind1),exp(siglevel_local_extr(ind1)),'color','g','marker','+');
    end
    plot(hax1,intersect_expfit(ind1),yield_cutlevel(ind1),'color','b','marker','+')
    str_leg{ind1} = [num2str(power_sorted(ind1)) ' mW; ' num2str(PeakIntensity(ind1)*1e-12) ' TW/cm^2; ' num2str(E_cutoff(ind1)) ' eV'];
    if 0
        figure;plot(x_expfit,y_expfit,'k')
        hold on;
        plot(intersect_expfit(ind1),FIT_tail(ind1).m(1)*exp(FIT_tail(ind1).m(2)*intersect_expfit(ind1)),'r+')
        showfit(FIT_tail(ind1))
        hold off;
    end
end

% xlim([0 E_max]);
ylim([min(abs(norm_minmax(ind_select2,1))) max(abs(norm_minmax(ind_select2,2)))]);
xlabel(hax1,'E_{kin} [eV]')
ylabel(hax1,'Counts/s (normalized)')
title(hax1,'Spectrum over laser power and cutoff fit')
% str_leg{1} = '';
legend(hndl1(ind_select2),str_leg(ind_select2));
if flag_save
    saveas(SPECTRA_DETAILS,[folder_to_save 'spectra_details.fig']);
end

if exist('params_for_analyze_FE_itx.m')==0
    if 1
        ind_select3 = 3:9;
    else
        ind_select3 = vec2ind((E_diff>=0).*ind2vec(ind_select2,Nspec));
    end
end

x_plot1 = PeakIntensity(ind_select2);
y_plot1 = E_cutoff(ind_select2) - 1*cikkfactor*WorkFunction_eV;
x_fit_FE = PeakIntensity(ind_fit_FE);
y_fit_FE = E_cutoff(ind_fit_FE) - 1*cikkfactor*WorkFunction_eV;
FIT_FE = linear_fit_to_data(x_fit_FE,y_fit_FE,E_cutoff_std(ind_fit_FE));
% slope_FE = FIT_FE.m(2)*10^11; %larger number for matlab
fit_errors_FE = lin_fit_error_calculator( transpose(x_fit_FE), transpose(y_fit_FE), FIT_FE.m(1), FIT_FE.m(2), transpose(E_cutoff_std_chi_sq(ind_fit_FE)) ); %calc error from results with chi square methode
E_cut_slope = FIT_FE.m(2)*1e-4*C_SI.qE; % [J/(W/m^2)]  
field_enh0 = sqrt(E_cut_slope*2*C_SI.mE*C_SI.eps0*C_SI.c*omega^2/C_SI.qE^2/10.007);
field_enh0_std = 1/2*sqrt(C_SI.mE*C_SI.eps0*C_SI.c*omega^2/C_SI.qE^2/10.007)*1/sqrt(E_cut_slope)*FIT_FE.stats.coeffs(2,2)*1e-4*C_SI.qE;
field_enh0_std_chi_sq = 1/2*sqrt(C_SI.mE*C_SI.eps0*C_SI.c*omega^2/C_SI.qE^2/10.007)*1/sqrt(E_cut_slope)*fit_errors_FE.error.slope*1e-4*C_SI.qE;

E_diff = (E_cutoff - FIT_FE.m(1) - cikkfactor*WorkFunction_eV);
E_diff2 = (E_cutoff - cikkfactor*WorkFunction_eV);
E_diff_std_chi_sq = E_cutoff_std_chi_sq + fit_errors_FE.error.axis_cr;
eplasmon = zeros(size(E_diff));
field_enh = zeros(size(E_diff));
eplasmon_std_chi_sq = zeros(size(E_diff));
eplasmon(ind_select3) = sqrt((E_diff(ind_select3)*eVtoJoule*16*pi^2*C_SI.mE*C_SI.c^2)/(10.007*C_SI.qE^2*lambda^2));
eplasmon_oldw(ind_select3) = sqrt((E_diff2(ind_select3)*eVtoJoule*16*pi^2*C_SI.mE*C_SI.c^2)/(10.007*C_SI.qE^2*lambda^2));
eplasmon_std_chi_sq(ind_select3) = 0.5*sqrt((16*pi^2*C_SI.mE*C_SI.c^2)/(10.007*C_SI.qE^2*lambda^2))./sqrt(E_diff(ind_select3)*eVtoJoule).*E_diff_std_chi_sq(ind_select3)*eVtoJoule;
field_enh(ind_select3) = eplasmon(ind_select3)./ElectricField(ind_select3);
% calculate the field enhancement without the cutoff scaling fit field_enh_oldw
field_enh_oldw(ind_select3) = eplasmon_oldw(ind_select3)./ElectricField(ind_select3);
field_enh_oldw(imag(field_enh_oldw) ~= 0) = 0; %set imagenary part to zero

field_enh_std_chi_sq(ind_select3) = eplasmon_std_chi_sq(ind_select3)./ElectricField(ind_select3);
if length(ind_select3)~=length(E_diff)
    disp('!!!');
    disp('Warning: Some values of the cutoff energy are too low for F.E. calculation!');
    disp('!!!');
end

disp(['peak intensity [W/cm^2] = ' num2str(map2rowvec(PeakIntensity(ind_select3)),'%1.3e\t')]);
disp(['cutoff energy [eV] = ' num2str(map2rowvec(E_cutoff(ind_select3)),'%1.3f\t')]);
disp(['laser field [V/m] = ' num2str(map2rowvec(ElectricField(ind_select3)),'%1.3e\t')]);
disp(['plasmon field [V/m] = ' num2str(map2rowvec(eplasmon(ind_select3)),'%1.3e\t')]);
disp(['field enhancement = ' num2str(map2rowvec(field_enh(ind_select3)),'%2.3f\t')]);
%% Plot cutoff scaling, FE and Keldysh parameters
CUTOFF = figure;
fontsize = 20;
linewidth = 1.5;
if exist('params_for_analyze_FE_itx.m')==0
    ind_fit_cutoff = [1 length(ind_select3)];
end

plot(x_plot1,y_plot1,'k+')
showfit(FIT_FE);
xlim([0 max(x_plot1)*1.1]);
ylim([min(0,min(y_plot1)) max(y_plot1)*1.1]);
title(['F.E. from slope = ' num2str(field_enh0)])
ylabel('Cutoff energy - 0.538\cdotW [eV]','FontSize',fontsize), xlabel('Peak intensity [W/cm^2]','FontSize',fontsize)
setfigP;
if flag_save
    saveas(CUTOFF,[folder_to_save 'cutoff_scaling.fig']);
end

if flag_plateau
    CUTOFF_error = figure;
    fontsize = 20;
    linewidth = 1.5;
    if exist('params_for_analyze_FE_itx.m')==0
        ind_fit_cutoff = [1 length(ind_select3)];
    end

    plot(x_plot1,y_plot1,'k+')
    errorbar(x_plot1,y_plot1,E_plateau_std_chi_sq(ind_select2), 'ko')
    showfit(FIT_FE);
    xlim([0 max(x_plot1)*1.1]);
    ylim([min(0,min(y_plot1)) max(y_plot1)*1.1]);
    title(['F.E. from slope = ' num2str(field_enh0)])
    ylabel('Cutoff energy - 0.538\cdotW [eV]','FontSize',fontsize), xlabel('Peak intensity [W/cm^2]','FontSize',fontsize)
    setfigP;
    if flag_save
        saveas(CUTOFF_error,[folder_to_save 'cutoff_scaling_err.fig']);
    end
end
    
FE_CUTOFF = figure('units','normalized','outerposition',[0 0 1 1]);
[hAx,hLine1,hLine2] = plotyy(PeakIntensity(ind_select3),field_enh(ind_select3),PeakIntensity(ind_select3),E_cutoff(ind_select3));
hLine1.LineStyle = 'none';
hLine2.LineStyle = 'none';
set(hLine1,'Marker','x','Linewidth',linewidth, 'MarkerSize', 14)
set(hLine2,'Marker','x','Linewidth',linewidth, 'MarkerSize', 14)
set(hAx,'Linewidth',linewidth,'FontSize',fontsize)

title(['F.E. from slope = ' num2str(field_enh0)])
xlabel('Peak intensity [W/cm^2]')
ylabel(hAx(1),'Field enhancement') % left y-axis 
ylabel(hAx(2),'Cutoff energy [eV]','FontSize',fontsize) % right y-axis
if flag_save
    saveas(FE_CUTOFF,[folder_to_save 'FE_and_cutoff_vs_intensity.fig']);
end
%%
keldysh0 = omega*sqrt(2*C_SI.mE*WorkFunction_eV*eVtoJoule)/C_SI.qE./field_enh0./ElectricField;
keldysh = omega*sqrt(2*C_SI.mE*WorkFunction_eV*eVtoJoule)/C_SI.qE./field_enh./ElectricField;
keldysh_oldw = omega*sqrt(2*C_SI.mE*WorkFunction_eV*eVtoJoule)/C_SI.qE./field_enh_oldw./ElectricField;
keldysh0_std_chi_sq = abs(-omega*sqrt(2*C_SI.mE*WorkFunction_eV*eVtoJoule)/C_SI.qE./(field_enh0^2)./ElectricField).*field_enh0_std_chi_sq;
keldysh_std_chi_sq = abs(-omega*sqrt(2*C_SI.mE*WorkFunction_eV*eVtoJoule)/C_SI.qE./(field_enh.^2)./ElectricField).*field_enh_std_chi_sq;


FE_KELDYSH0 = figure;
[ax1,hl1,hl2] = plotxxP(PeakIntensity(ind_select3),keldysh0(ind_select3),field_enh(ind_select3),[{'Incident Peak Intensity [W/cm^2]'},{'Local Keldysh parameter'}],{'Field Enhancement'});
set(ax1,'Linewidth',linewidth,'FontSize',fontsize);
legend(['F.E. from slope = ' num2str(roundP(field_enh0,1))  char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
title('Keldysh \gamma from slope');
if flag_save
    saveas(FE_KELDYSH0,[folder_to_save 'FE_vs_Keldysh0.fig']);
end

FE_KELDYSH = figure;
[ax1,hl1,hl2] = plotxxP(PeakIntensity(ind_select3),keldysh(ind_select3),field_enh(ind_select3),[{'Incident Peak Intensity [W/cm^2]'},{'Local Keldysh parameter'}],{'Field Enhancement'});
set(ax1,'Linewidth',linewidth,'FontSize',fontsize);
legend(['F.E. from slope = ' num2str(roundP(field_enh0,1))  char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
title('Keldysh parameter calculated separately for every intensity');
if flag_save
    saveas(FE_KELDYSH,[folder_to_save 'FE_vs_Keldysh.fig']);
end

FE_KELDYSH0_err = figure;
% [ax1,hl1,hl2] = plotxxP(PeakIntensity(ind_select3),keldysh(ind_select3),field_enh(ind_select3),[{'Incident Peak Intensity [W/cm^2]'},{'Local Keldysh parameter'}],{'Field Enhancement'});
subplot(1,2,1)
errorbar(keldysh0(ind_select3),field_enh(ind_select3), field_enh_std_chi_sq(ind_select3),  'ko', 'Linewidth',linewidth)
set(gca,'Linewidth',linewidth,'FontSize',fontsize);
legend(['F.E. from slope = ' num2str(roundP(field_enh0,1)) char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
xlabel('keldysh')
ylabel('Field enhancement')
subplot(1,2,2)
errorbar(keldysh0(ind_select3),field_enh(ind_select3), keldysh0_std_chi_sq(ind_select3),'horizontal', 'ko', 'Linewidth',linewidth)
set(gca,'Linewidth',linewidth,'FontSize',fontsize);
% legend(['F.E. from slope = ' num2str(roundP(field_enh0,1)) char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
xlabel('Keldysh \gamma')
ylabel('Field enhancement')
title('Keldysh \gamma from slope');
if flag_save
    saveas(FE_KELDYSH0_err,[folder_to_save 'FE_vs_Keldysh0_err.fig']);
end

FE_KELDYSH_err = figure;
% [ax1,hl1,hl2] = plotxxP(PeakIntensity(ind_select3),keldysh(ind_select3),field_enh(ind_select3),[{'Incident Peak Intensity [W/cm^2]'},{'Local Keldysh parameter'}],{'Field Enhancement'});
subplot(1,2,1)
errorbar(keldysh(ind_select3),field_enh(ind_select3), field_enh_std_chi_sq(ind_select3),  'ko', 'Linewidth',linewidth)
set(gca,'Linewidth',linewidth,'FontSize',fontsize);
legend(['F.E. from slope = ' num2str(roundP(field_enh0,1)) char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
xlabel('keldysh')
ylabel('Field enhancement')
subplot(1,2,2)
errorbar(keldysh(ind_select3),field_enh(ind_select3), keldysh_std_chi_sq(ind_select3),'horizontal', 'ko', 'Linewidth',linewidth)
set(gca,'Linewidth',linewidth,'FontSize',fontsize);
% legend(['F.E. from slope = ' num2str(roundP(field_enh0,1)) char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
xlabel('Keldysh parameter')
ylabel('Field enhancement')
title('Keldysh parameter calculated separately for every intensity');
if flag_save
    saveas(FE_KELDYSH_err,[folder_to_save 'FE_vs_Keldysh_err.fig']);
end


FE_KELDYSH_oldw = figure;
[ax1,hl1,hl2] = plotxxP(PeakIntensity(ind_select3),keldysh0(ind_select3),field_enh_oldw(ind_select3),[{'Incident Peak Intensity [W/cm^2]'},{'Local Keldysh parameter'}],{'Field Enhancement'});
set(ax1,'Linewidth',linewidth,'FontSize',fontsize);
% legend(['F.E. from slope = ' num2str(roundP(field_enh0,1))  char(177) num2str(roundP(field_enh0_std_chi_sq, 3))])
legend(['F.E. mean = ' num2str(roundP(mean(field_enh_oldw(4:9)),1)) ' ,F.E. std = ' num2str(roundP(std(field_enh_oldw(4:9)),1)) ' (points:4-9)'])
title('F.E. calculated without cutoff scaling fit');
if flag_save
    saveas(FE_KELDYSH_oldw,[folder_to_save 'FE_vs_Keldysh_oldw.fig']);
end

if flag_plateau==0
    disp('!!! For error calculation of the tail fit, error calculation for the exponential fit needs to be implemented!');
end

%% variance of cutoff scaling fit

% determine the the error for "a" and "b" (a*x+b)
% assumption: x values are precise(no error)

% fit_errors = lin_fit_error_calculator( x_linfit1, y_linfit1, FIT_sp(ind1,3).m(1), FIT_sp(ind1,3).m(2))

% x_data = x_fit_FE;
% y_data = y_fit_FE;
% 
% a = FIT_FE.m(1);
% b = FIT_FE.m(2);
% % 
% % 
% %     sigma_i = ones(length(x_data),1);
% sigma_i = E_cutoff_std_chi_sq(ind_fit_FE);
%     s = sum(sigma_i.^2);
%     s_x = sum(x_data./(sigma_i.^2));
%     s_xx = sum((x_data.^2)./(sigma_i.^2));
%     delta = s*s_xx-(s_x)^2;
% %     chi_square = sum(((y_data-(a+x_data*b)).^2)./sigma_i);
% %     sigma_square = chi_square/(length(x_data)-2);
%     sigma_a_square = (s_xx/delta);
%     sigma_b_square = (s/delta);


%% Save to .mat file
struct_data.raw.data_raw = data_raw;
struct_data.raw.nrg = nrg;
struct_data.raw.alpha = alpha;
struct_data.raw.power = power;

struct_data.int_calib.focus_x_um = focus_x_um;
struct_data.int_calib.focus_y_um = focus_y_um;
struct_data.int_calib.RepetitionRate = RepetitionRate;
struct_data.int_calib.PulseLength_fs = PulseLength_fs;
struct_data.int_calib.lambda_nm = lambda*1e9;
struct_data.int_calib.WorkFunction_eV = WorkFunction_eV;

struct_data.params.flag_plateau = flag_plateau;
struct_data.params.flag_E0_from_fit = flag_E0_from_fit;
struct_data.params.peak_part = peak_part;
struct_data.params.E_fit0 = E_fit0;
struct_data.params.E_fit1 = E_fit1;
struct_data.params.E_fit2 = E_fit2;
struct_data.params.E_fit3 = E_fit3;
struct_data.params.ind_select0 = ind_select0;
struct_data.params.ind_select1 = ind_select1;
struct_data.params.ind_select2 = ind_select2;
struct_data.params.ind_select3 = ind_select3;
struct_data.params.ind_fit_int = ind_fit_int;
struct_data.params.ind_fit_cutoff = ind_fit_cutoff;

struct_data.process.Edata_plot = Edata_plot;
struct_data.process.Edata_n_sigma = Edata_n_sigma;
struct_data.process.power_sorted = power_sorted;
struct_data.process.E_maxyield0 = E_maxyield0;
struct_data.process.E0 = E0;
struct_data.process.E0_std = E0_std;
struct_data.process.E_plateau = E_plateau;
struct_data.process.E_plateau_std = E_plateau_std;
struct_data.process.E_high = E_high;
struct_data.process.E_high_std = E_high_std;
struct_data.process.FIT_sp = FIT_sp;

struct_data.result.PeakIntensity = PeakIntensity(ind_select3);
struct_data.result.keldysh = keldysh(ind_select3);
struct_data.result.totalCounts = totalCounts(ind_select3);
struct_data.result.E_cutoff = E_cutoff(ind_select3);
struct_data.result.E_cutoff_std = E_cutoff_std(ind_select3);
struct_data.result.ElectricField = ElectricField(ind_select3);
struct_data.result.eplasmon = eplasmon(ind_select3);
struct_data.result.field_enh = field_enh(ind_select3);
struct_data.result.field_enh0 = field_enh0;
struct_data.result.field_enh0_std = field_enh0_std;
if flag_save
    save([folder_to_save 'data.mat'],'struct_data');
    fig2png(folder_to_save);
end
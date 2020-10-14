purge;
%% Settings for the run
save_only_essentials = 0;
num_retrievals = 1;

% Choose how plotting interval is determined
if 0 % Plot every X seconds
   plotting_criterion = 'time';
   plotting_interval = 60*30; % seconds
elseif 1 % Plot every Nth generation
   plotting_criterion = 'generation';
   plotting_interval = 10;
else % Plot after every Nth evaluation
   plotting_criterion = 'evaluation';
   plotting_interval = 1000;
end
%% Trace and DE Parameters (Default settings, controlled individually for each experiment)
E_length = 512; % this parameter is proportional to the time window considered if 'dt_in_cycles' parameter is kept constant.
dt_in_cycles = 0.2;
tau_int_selection = 'all';
tau_lims = [-1,1]*120;

time_limit = 1*60*60;           %time limit (in seconds)
pop_size = 100;                  %population size
obj_fun_eval_limit = 2e2*pop_size;       %objective function evaluation limit
CR = 0.5;                       %crossover probability
num_points_start = E_length/2;	%starting number of points
num_points_max = E_length;             %maximum number of points - VERY IMPORTANT PARAMETER!!! - have to be kept sufficiently high or 0
max_sampling_in_cycles = 1/5; % only used if num_points_max is set to 0
resolution_rate = 50;           %resolution increases every Nth generation
convergence_condition = 0e-5;   %detects convergence based on standard deviation of the latest generation

random_sample_points = 1;
shift_peak_to_zero_delay = 1;
mirror_bands = 0;               %choose +-1 to mirror the subtrace, this ensures symmetry with respect to delay. 0 disables mirroring.
survival_of_the_fittest = 1;
gaussian_filtering = 0;
% Adds normally distributed noise to measured trace. The noise factor is
% the standard deviation of the noise distribution, with respect to peak
% intensity of the FROG trace that is scaled to unity.
noise_factor_trace = 0*5e-2;
%% Experimental data: Load & process trace
input_read_IFROG.path_dir = 'FROG4\';
input_read_IFROG.save_fig = 0;  % save generated figures and data struct at the end
flag_load_seed = 1;
pulse_opt.test = 0;
flag_load_spectrum = 1;         % load spectrum for plotting it overlaid with reconstructed spectrum
input_DE_IFROG.clean_spectrum = 1;  % remove spectral components outside the measured laser bandwidth at every iteration
input_DE_IFROG.vec_GDD = [0 25 50 75 100];
input_DE_IFROG.vec_TOD = [-40 -20 0 20 40];
band_selection = [1 1 0];
weight = [1 1 0];

input_read_IFROG.delay_select = [0 0];
input_read_IFROG.ignore_delays_start_end = [0 0];
input_read_IFROG.use_wavelengths = [330 525]; % [nm]
input_read_IFROG.artefact_removal_iterations = 0;
input_read_IFROG.artefact_tol = 1.5;
input_read_IFROG.remove_background = 1;
% method_for_finding_trace_center_point = 'harmonic center of gravity';
input_read_IFROG.method_for_finding_trace_center_point = 'use margins';
input_read_IFROG.remove_fundamental_bg = 0;
input_read_IFROG.clean_periphery_opt = 0;
input_DE_IFROG.save_fig = input_read_IFROG.save_fig;
input_read_IFROG.plot_log = 0;

file_spectrum = 'Venteon_spectrum_norm.txt';
% file_spectrum = 'U:\People\S.Peter\projects\IFROG\Venteon_spectrum_from_datasheet.txt';
% file_spectrum = 'U:\Measurement_Data\OAC_IFROG\2020_03_10_ACER\Venteon_spectrum_noheader.txt';
seed_sp_amp = ['E_vs_omega.txt'];
if flag_load_seed
    input_DE_IFROG.seed_sp_amp = seed_sp_amp;
else
    input_DE_IFROG.seed_sp_amp = [];
end

if pulse_opt.test
    temp = load(seed_sp_amp);
    pulse_opt.t_L = temp(:,1);
    pulse_opt.Emag = abs(temp(:,2) + 1i*temp(:,3));
    pulse_opt.Ephase = unwrap(angle(temp(:,2) + 1i*temp(:,3)));
end

if flag_load_spectrum
    file_spectrum = 'U:\Measurement_Data\OAC_IFROG\2020_10_06_ACER\Venteon_spectrum.txt';
    data_sp = read_spectrum(file_spectrum);
    input_DE_IFROG.fund_spectrum = condition_spectrum(data_sp);
    l0 = sum(input_DE_IFROG.fund_spectrum(:,1).*input_DE_IFROG.fund_spectrum(:,2))/sum(input_DE_IFROG.fund_spectrum(:,2))*1e-9;
else
    input_DE_IFROG.fund_spectrum(:,1) = 2*pi*map2colvec(0.360:0.0005:0.390);
    input_DE_IFROG.fund_spectrum(:,2) = exp(-(input_DE_IFROG.fund_spectrum(:,1)-2*pi*0.375).^2/0.02^2);
end

% Determine supergaussian filter parameters based on independently
% measured fundamental spectrum
if input_DE_IFROG.clean_spectrum
    fund_spectrum = load(file_spectrum);
    fund_spectrum(:,2) = (fund_spectrum(:,2)-mean(fund_spectrum(1:50,2)))/max(fund_spectrum(:,2)-mean(fund_spectrum(1:50,2)));
    fund_spectrum((fund_spectrum(:,2)<0.003),2) = 0;
    sp_om = SpectrumConvert(fund_spectrum);
    omega0 = peak_props([sp_om(:,1) sp_om(:,2).^2]);
    index_nonzero = logical(sp_om(:,2)>0);
    omega0_SG = (min(sp_om(index_nonzero,1)) + max(sp_om(index_nonzero,1)))/2;
    order_SG = 8;
    sigma_SG = fit_SG_sigma(sp_om(:,1),index_nonzero,order_SG,omega0_SG); % determine width of supergaussian filter
    input_DE_IFROG.SG_params = [order_SG sigma_SG omega0_SG]; % parameters for supergaussian filter: [order, sigma, omega0]
    if 1
        hfig_sp_filter = figure;hold on;
        plot(sp_om(:,1),sp_om(:,2)/max(sp_om(:,2)),'k');
        plot(sp_om(:,1),sp_om(:,2)/max(sp_om(:,2)).*exp(-((sp_om(:,1)-omega0_SG).^2/sigma_SG^2).^order_SG),'k--');
        plot(sp_om(:,1),exp(-((sp_om(:,1)-omega0_SG).^2/sigma_SG^2).^order_SG),'r');
        xlabel('\omega [rad/fs]')
        legend('reference spectrum','filter applied','filter function')
    end
end

% Spectral weighting
BBO = 0;
pm_min_level = 0.0;
input_read_IFROG.Protected_Ag_mirror = 0;
input_read_IFROG.BG23 = 0;
input_read_IFROG.BG39 = 0;
input_read_IFROG.UV_focusing_lens = 0;
input_read_IFROG.UV_enhanced_Al_mirror = 0;
input_read_IFROG.HR270 = 0;
input_read_IFROG.IF270 = 0;
input_read_IFROG.IFROG_order = 2;

% SHG DE-IFROG EXPERIMENTS
files1 = dir([input_read_IFROG.path_dir '*.txt']);
if exist([input_read_IFROG.path_dir 'FROG_dot.txt'])~=2
    for ind1 = 1:length(files1)
        comma2dot([input_read_IFROG.path_dir files1(ind1).name]);
    end
end

if ~strcmp(input_read_IFROG.path_dir(end),'\')
    input_read_IFROG.path_dir = [input_read_IFROG.path_dir '\'];
end
ind_path = 0;
while exist([input_read_IFROG.path_dir 'retrieval' num2str(ind_path) '\'])==7
   ind_path = ind_path+1;
end
% ind_path = 4;
input_read_IFROG.savepath1 = [input_read_IFROG.path_dir 'retrieval' num2str(ind_path) '\'];
input_DE_IFROG.savepath1 = input_read_IFROG.savepath1;
if input_read_IFROG.save_fig
    mkdir(input_read_IFROG.savepath1);
end


fiber_pick = 'FC-UV600';
grating_pick = 600;
l_center_pick = 415;

%%%% Calibration options
if 0 % Spectrometer + fiber calibration
 input_read_IFROG.calibration_opt = struct(...
    'fiber',fiber_pick,...
    'l_center',l_center_pick,...
    'grating',grating_pick);
 %          fiber_selection = {'FC-UV600','FC-UV100','QP400'};
 %          grating_selection = [600, 1200];
 %          l_center_selection = [270, 415, 400, 385];

%          calibration_opt = struct(...
%             'fiber',fiber_selection(fiber_pick),...
%             'l_center',l_center_selection(l_center_pick),...
%             'grating',grating_selection(grating_pick));
else
 input_read_IFROG.calibration_opt = 0;
end
if BBO % Phase matching
 crystal_length = 20e-6;
 crystal_angle  = 29.2;
 input_read_IFROG.phasematching_opt = struct(...
    'crystal_length',crystal_length,...
    'crystal_angle',crystal_angle,...
    'min_level',pm_min_level);
else
 input_read_IFROG.phasematching_opt = 0;
end

if exist([input_read_IFROG.path_dir 'S.mat'])==2
    temp = load([input_read_IFROG.path_dir 'S.mat']); S = temp.S;
    input_DE_IFROG.name = S.name;
else
    S = Read_IFROG_data(input_read_IFROG);
    save([input_read_IFROG.path_dir 'S.mat'],'S');
end

% clear fiber_pick fiber_pick...
%  grating_pick grating_pick...
%  l_center_pick l_center_pick...
%  path_dir plot_log save_fig;

% Use data set
T_f = S.T_f;
tau = S.tau;
f = S.f;
if ~flag_load_spectrum
    l0 = S(1).l0;
end
%% Synthetic data: Use test pulse
if 0
   T_given = 0; tau_given = 0; f_given = 0;
   l0 = 800e-9;
   pulse_width = 40e-15;
   
   GDD = 1*27.796*3.0; %fs^2/rad CaF_2 = 27.796fs^2/mm, 3.04mm
   TOD = 0; %fs^3/rad^2 300
   FOD = 0; %fs^4/rad^3 1400
   SPM = 1/4*pi; % rad?
   phase_noise = 0; % 0 --> no random factor
   
   use_predefined_test_phase = 0;
   
   % Satellites are given as an array of size (N,3), with each row defined as:
   % [amplitude (x (times pulse amplitude),
   %  width (x times pulse width),
   %  position (x times pulse width)]
   % satellites = [0.5,1.0,-3;...
   %               0.3,0.5,+2];
   % satellites = [0.4,0.9,+1.7];
   satellites = [0.6,0.9,+2;
       0.3,0.8,-2];
   % satellites = 0;
   
   pulse_opt = struct(...
      'pulse_width',pulse_width,...
      'GDD',GDD,...
      'TOD',TOD,...
      'FOD',FOD,...
      'SPM',SPM,...
      'phase_noise',phase_noise,...
      'satellites',satellites,...
      'use_predefined_test_phase',use_predefined_test_phase);
   
   % Trace parameters
   IFROG_order = 2;
   band_selection = [1 1 0];
   weight = [1 1 0];
%    tau_lims = [-1,1]*180; %176
   c = 299792458;
   dt = (l0/c) * 1.0;
   
   name = '00000tenpercentnoisetrace';
else
   T_given = T_f;
   f_given = f;
   tau_given = tau;
end
%% Reference spectrum
if 0
   %    filename_ref = '20170214_SPIDER_reference_spectrum_for_DE-IFROG_CORRECT.mat';
   filename_ref = '20170405_reference_spectrum_for_DE-IFROG_AFTER_IFROG_FOCUS_WITH_PAPER.mat';
   S_ref_spectrum = load([Path_MBI_projects ...
      '2016_IFROG_genetic\Data\20170214\SPIDER\' filename_ref]);
   %    S_ref_spectrum = load('Z:\Projektit\2017_MISC\2017-02-10_SPIDER_flat_plus_1m_air.mat');
   S_ref_spectrum = S_ref_spectrum.S_ref;
elseif 0 % Load Konstanz Fundamental spectrum      
   A = dlmread([Path_MBI_projects ...
      '2017_THIFROG_Konstanz\Data\Marginal_correction\' ...
      'Konstanz_marginal_correction_data.txt']);
   l_fund = Make_row(A(:,1)*1e-9);
   I_l_fund = Make_row(A(:,2));
   I_l_fund_edit = I_l_fund;
   I_l_fund_edit(I_l_fund_edit<0)=0;
   I_l_fund_edit(1:Find_index(l_fund,850*1e-9))=0;
   I_l_fund_edit = Make_row(smooth(I_l_fund_edit,3));
%    figure
%    plot(l_fund,I_l_fund,l_fund,I_l_fund_edit)
   
   I_l_fund = Normalize(I_l_fund_edit);
   I_f_fund = flip(Normalize(I_l_fund.*l_fund.^2));
   c = 299792458;
   f_fund = flip(c./l_fund);
   S_ref_spectrum = struct;
   S_ref_spectrum.f = Make_row(f_fund);
   S_ref_spectrum.I_f = Make_row(I_f_fund);
   S_ref_spectrum.phi_f = zeros(size(S_ref_spectrum.f));
   S_ref_spectrum.l = Make_row(l_fund);
   S_ref_spectrum.I_l = Make_row(I_l_fund);
   S_ref_spectrum.phi_l = zeros(size(S_ref_spectrum.l));
else
   S_ref_spectrum = 0;
end
%% DE
% Compute delay step
c = 299792458;
dt = (l0/c) * dt_in_cycles;

% If maximum number of points was not given, estimate it
if num_points_max == 0
   num_points_max = round(dt * E_length/2/(l0/c)*max_sampling_in_cycles);
end

if 1 % Parameters into structs
   T_opt = struct(...
      'IFROG_order',input_read_IFROG.IFROG_order,...
      'band_selection',band_selection,...
      'weight',weight,...      
      'E_length', E_length,...
      'tau_int_selection',tau_int_selection,...
      'tau_lims',tau_lims,...
      'dt',dt,...
      'l0',l0,...
      'T_given',T_given,...
      'tau_given',tau_given,...
      'f_given',f_given,...
      'noise_factor_trace',0); %noise_factor_trace);
   
   DE_opt = struct(...
      'CR',CR,...
      'pop_size',pop_size,...
      'time_limit',time_limit,...
      'obj_fun_eval_limit',obj_fun_eval_limit,...
      'num_points_start',num_points_start,...
      'num_points_max',num_points_max,...
      'resolution_rate',resolution_rate,...
      'plotting_interval',plotting_interval,...
      'plotting_criterion',plotting_criterion,...
      'convergence_condition',convergence_condition,...
      'survival_of_the_fittest',survival_of_the_fittest,...
      'mirror_bands',mirror_bands,...
      'shift_peak_to_zero_delay',shift_peak_to_zero_delay,...
      'gaussian_filtering',gaussian_filtering,...
      'random_sample_points',random_sample_points);
end


input_DE_IFROG.save_only_essentials = 0;
% num_retrievals = 1;
% C_R = cell(num_retrievals,1);
% for cc = 1:num_retrievals
% close all
Disp_char_line('*',3)
Disp_title(['Starting retrieval'])
Disp_char_line(' ')
Disp_char_line('*',3)
%     cc
T_opt.noise_factor_trace = 1E-3;

C_R = DE_IFROG(T_opt, DE_opt, pulse_opt, S_ref_spectrum, input_DE_IFROG); % Main function

if input_DE_IFROG.save_fig
    saveas(hfig_sp_filter,[input_DE_IFROG.savepath1 'spectral_filter.fig'])
    save([input_DE_IFROG.savepath1 'C_R.mat'],'C_R');
end
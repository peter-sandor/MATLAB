% *******************************************
% Weinacht Group Strong Field Ionization Code
% *******************************************
% *******************************************
% This code solves for the time-dependent state
% populations & photoelectron spectrum for a
% femtosecond laser pulse driving an inverted A system:
% 
% /////// Continua ///////
%   ^                ^
%   |(n12 photons)   |(n34 photons)
%   |                |
% ---- |1>  <--->  ---- |3>
%             ^
%             |   (n01 photons)
%             |
%           ----- |0>
% 
% *******************************************
% Start from a single ground state, have two intermediate bound states,
% which are coupled and two separate continua
% *******************************************
function output = SFI_Inverted_A_for_paper(input)
% SECTION I: System and Laser Parameters (SI units, all times in ps)
% A: Constants
params.HBAR = 1.054572E-34; %Reduced Planck's constant,[J*s]
params.ALPHA = 1*1/(2*params.HBAR*2)*1e-12; %Coupling parameter for continuum ladder
params.EMAX = 2*1.6e-19; %Maximum photoelectron energy [J] (i.e., cutoff)
params.M = 400; %Number of states in each continuum ladder
NUMBER_OF_BOUND_STATES = 3; %Must match number of bound states 

% B: Laser
params.wField = 2*pi*input.nu_field; %Central frequency of laser pulse, [radian/ps]
FWHM = input.FWHM; %FWHM of laser pulse, [ps]
params.Io = input.intensity; %Intensity of laser pulse, [W/m^2]
params.I_REF = 1.278e17; %Reference intensity to which peak Rabi frequencies and Stark shifts (defined below) correspond, [W/m^2] 
% (Add e.g. pulse shape parameters here)

% C: Levels
params.w0 = 0; %Peak AC Stark shift of level |0>, [radian/ps]
params.w1 = 1*2*pi*180; % 184 THz %Peak AC Stark shift of level |1>
params.w3 = 1*2*pi*180; %Peak AC Start shift of level |3>
params.wC = 1*2*pi*180; %Peak AC Stark shift of Continuum
% D: Transitions
%   Notation is (QUANTITY)(Lower Level)(Upper Level)
%   w = transition frequency ; n = photon number;
%   rabi = peak amplitude of transition matrix element
%   e.g. rabi01 is rabi frequency for |0>->|1> transition

params.n01 = 5;
params.n1C1 = 2;
params.n0C1 = params.n01 + params.n1C1;
params.n03 = 5;
params.n3C3 = 2;
params.n0C3 = params.n03 + params.n3C3;

params.w01 = 2*pi*1860; % 1800 THz (1860)
params.w03 = 2*pi*1920; % 1825 THz (1920)
params.w0C3 = 2*pi*2495 + params.EMAX/params.HBAR/2*1E-12; % 2478 THz (10.26 eV) for D1 (See paper for info on the params.EMAX term)
params.w0C1 = 2*pi*2340 + params.EMAX/params.HBAR/2*1E-12; % 2340 THz (9.69 eV) for D0 (See paper for info on the params.EMAX term)
% The energy axis calibration depends on the EMAX parameter. Also, the way
% the calculation framework is set up, it includes an energy shift of
% EMAX/2 added to what the peak location should be. To correct for that,
% either the energy axis can be adjusted or the ionization potential raised by
% EMAX/2.
params.w1C1 = params.w0C1-params.w01;
params.w3C3 = params.w0C3-params.w03;

%Next lines are the Rabi Frequencies between states |0>,|1>,|3>,|C1> and |C3>
params.rabi0C1 = 0*2*pi*3;
params.rabi01 = 1*2*pi*1;
params.rabi1C1 = 2*2*pi*10;

params.rabi0C3 = 0*2*pi*3;
params.rabi03 = 1*2*pi*10;
params.rabi3C3 = 2*2*pi*10;

params.coupling13 = 1/2*2*pi*15; % Nonadiabatic coupling between 1 and 3
% % II: Laser Pulse Shaping

% A: Define time axis from -1000fs to +1000fs
step=0.001;
lim = 1; 
times=-lim:step:lim;

% B: Define an unshaped Gaussian pulse with specified FWHM
tau = 2*FWHM/2.35;
e_t=exp(-times.^2/tau^2);

% C: Fourier transform the pulse
% e_f=fft(fftshift(e_t));
e_f=fft(e_t);

% D: Apply frequency-domain SHAPING, if desired
f=fftshift(-1/(2*step):1/(2*lim):1/(2*step)); % Define frequency axis, centered about f = 0
phase=input.GDD/2*(2*pi*f).^2+input.TOD/6*(2*pi*f).^3+pi/2*sign(f+input.pistep); % Define a spectral phase profile; here, a pi phase jump at the central frequency
e_f_shaped=e_f.*exp(1i*phase); %Apply spectral phase profile

% E: Transform shaped pulse back to time domain
e_t_shaped=ifft(e_f_shaped);
timePulse = e_t_shaped; %Complex Electric Field vs. Time


% *******************************************

% III: Hamiltonian (See Paper for Discussion)

% This function represents the RHS of the TDSE written as:
%    Dt(s) = -i/hbar . H . s
% but is loosely named "Hamiltonian" in this code

function sp=Hamiltonian(t,s)

N = 2*params.M + NUMBER_OF_BOUND_STATES; %Number of rows in Hamiltonian
ind2=NUMBER_OF_BOUND_STATES+params.M; %Index where second bound state coupled to second continuum starts
% Preliminary: Interpolation of the Pulse E Field
ii=find(times>t,1);
dt=times(ii)-times(ii-1);
g=(abs(timePulse(ii))^2-abs(timePulse(ii-1))^2)/dt*(t-times(ii-1))+abs(timePulse(ii-1))^2;%Interpolating the |E|^2 values
phi=(angle(timePulse(ii))-angle(timePulse(ii-1)))/dt*(t-times(ii-1))+angle(timePulse(ii-1));%Inperpolating the ARG(E) values

% A: First, set continuum ladder couplings in EVERY row
%states S(3)->S(402) are one continuum
%states S(403)->S(803) are second continuum
sp=-1i*params.wC*g*s*params.Io/params.I_REF; %Each level gets a ponderomotive shift...
sp(1:(ind2-2))=sp(1:(ind2-2))-1i*params.EMAX*params.ALPHA*s(2:ind2-1);
sp(2:ind2-1)=sp(2:ind2-1)-1i*params.EMAX*params.ALPHA*s(1:(ind2-2));
sp(ind2+1:N-1)=sp(ind2+1:N-1)-1i*params.EMAX*params.ALPHA*s(ind2+2:N);
sp(ind2+2:N)=sp(ind2+2:N)-1i*params.EMAX*params.ALPHA*s(ind2+1:N-1);
%sp(1:(N-1))=sp(1:(N-1))-1i*params.EMAX*params.ALPHA*s(2:N); %each guy in the ladder is coupled to the one above and below it...
%sp(2:N)=sp(2:N)-1i*params.EMAX*params.ALPHA*s(1:(N-1)); 

%First 402 states are first two bound states and then 400 continuum ladder
%states. Next 401 states are the second bound state \3> above and the other
%400 continuum ladder states
% sp(Index_For_Cont_2)=-1i*wC*g*s(Index_For_Cont_2)*Io/I_REF;
% sp(Index_For_Cont_2)=sp(Index_For_Cont_2)-1i*params.EMAX*params.params.ALPHA*s(Index_For_Cont_2+1);
% sp(Index_For_Cont_2)=sp(Index_For_Cont_2)-1i*params.EMAX*params.ALPHA*s(NUMBER_OF_BOUND_STATES+2); %Need to couple to S(5)
% B: Overwrite bound state rows, leaving ladder couplings everywhere else
sp(1:3) = -1i*[s(1)*params.w0*g*params.Io/params.I_REF + s(2)*params.rabi01*(params.Io*g/params.I_REF)^(params.n01/2)*exp(-1i*(params.w01 - params.n01*params.wField)*t)*exp(-params.n01*1i*phi) + s(ind2)*params.rabi03*(params.Io*g/params.I_REF)^(params.n03/2)*exp(-1i*(params.w03 - params.n03*params.wField)*t)*exp(-params.n03*1i*phi) + params.rabi0C3*(params.Io*g/params.I_REF)^(params.n0C3/2)*exp(-1i*(params.w0C3 - params.n0C3*params.wField)*t)*s(ind2+1)*exp(-params.n0C3*1i*phi) + params.rabi0C1*(params.Io*g/params.I_REF)^(params.n0C1/2)*exp(-1i*(params.w0C1 - params.n0C1*params.wField)*t)*s(3)*exp(-params.n0C1*1i*phi),...
s(2)*params.w1*g*params.Io/params.I_REF + s(1)*params.rabi01*(params.Io*g/params.I_REF)^(params.n01/2)*exp(1i*(params.w01 - params.n01*params.wField)*t)*exp(params.n01*1i*phi) + params.rabi1C1*(params.Io*g/params.I_REF)^(params.n1C1/2)*exp(-1i*(params.w1C1 - params.n1C1*params.wField)*t)*s(3)*exp(-params.n1C1*1i*phi)+1i*params.coupling13*s(ind2),...
s(3)*params.wC*g*params.Io/params.I_REF + params.rabi1C1*(params.Io*g/params.I_REF)^(params.n1C1/2)*exp(1i*(params.w1C1 - params.n1C1*params.wField)*t)*s(2)*exp(params.n1C1*1i*phi)  + params.rabi0C1*(params.Io*g/params.I_REF)^(params.n0C1/2)*exp(+1i*(params.w0C1 - params.n0C1*params.wField)*t)*s(1)*exp(params.n0C1*1i*phi) + s(4)*params.EMAX*params.ALPHA,...
];
sp(ind2:ind2+1) = -1i*[s(ind2)*params.w3*g*params.Io/params.I_REF + s(1)*params.rabi03*(params.Io*g/params.I_REF)^(params.n03/2)*exp(1i*(params.w03 - params.n03*params.wField)*t)*exp(params.n03*1i*phi) + params.rabi3C3*(params.Io*g/params.I_REF)^(params.n3C3/2)*exp(-1i*(params.w3C3 - params.n3C3*params.wField)*t)*s(ind2+1)*exp(-params.n3C3*1i*phi)-1i*params.coupling13*s(ind2),...
s(ind2+1)*params.wC*g*params.Io/params.I_REF + params.rabi3C3*(params.Io*g/params.I_REF)^(params.n3C3/2)*exp(1i*(params.w3C3 - params.n3C3*params.wField)*t)*s(ind2)*exp(params.n3C3*1i*phi) + params.rabi0C3*(params.Io*g/params.I_REF)^(params.n0C3/2)*exp(+1i*(params.w0C3 - params.n0C3*params.wField)*t)*s(1)*exp(params.n0C3*1i*phi) + s(ind2+2)*params.EMAX*params.ALPHA,...
];
end
% *******************************************

% IV: Numerical Integration

% A: Set Desired Precision for Integration
PREC_ABS = 5E-6; %Absolute Precision
PREC_REL = 5E-6; %Relative Precision

% B: Set Initial State Vector and Call Integrator

systemSize = NUMBER_OF_BOUND_STATES + 2*params.M;
s0=zeros([systemSize 1]); %Create initial state vector
s0(1)=1; %All population initially in ground state

T_LIMIT = 0.2; %Limits of integration are +/-T_LIMIT
% tic
options = odeset('RelTol',PREC_ABS,'AbsTol',PREC_REL); %Options for integrator
[t,s] = ode23(@Hamiltonian,[-T_LIMIT T_LIMIT],s0,options); %Call integrator

% *******************************************

% V: Photoelectron Spectrum

% A: Select the Continuum Coefficients at Some Time
DESIRED_TIME = 0.1; %(Time at which to compute photoelectron spectrum)
s_ion1=s(:,NUMBER_OF_BOUND_STATES:ind2-1); %Extract the continuum state coefficients from the state vector
s_ion2=s(:,ind2+1:N);
CLOSEST_TIME = find(t > DESIRED_TIME,1); %(Picking the timestep closest to desired time)
ionPlot1 = s_ion1(CLOSEST_TIME,1:params.M); %The array of first continuum coefficients at CLOSEST_TIME
ionPlot2 = s_ion2(CLOSEST_TIME,1:params.M);

% B: Generate Photoelectron Spectrum Using Legendre Polynomials
eMaxReduced = params.EMAX/(1.602e-19); %Write params.EMAX in eV
resolution = 0.01; % [eV]
energies=0:resolution:eMaxReduced; % Create energy axis for the calculation

% For an explanation of the lines below, see the manuscript and
% references therein:

TEMP = size(ionPlot1);
MAX_L = TEMP(2);
runningSum1 = 0*double(energies);
runningSum2 = 0*double(energies);

for l=1:1:MAX_L
    legendAll = legendre(l-1,2*energies/eMaxReduced - 1);
    polynomial = legendAll(1,:);
    runningSum1 = ionPlot1(l)*sqrt((2*l-1)/eMaxReduced)*polynomial + runningSum1;
    runningSum2 = ionPlot2(l)*sqrt((2*l-1)/eMaxReduced)*polynomial + runningSum2;
end

PHOTOELECTRON_SPECTRUM = abs(runningSum1).^2+abs(runningSum2).^2; %The photoelectron spectrum

output.rsum(:,1)=abs(runningSum1).^2;
output.rsum(:,2)=abs(runningSum2).^2;
output.PES=PHOTOELECTRON_SPECTRUM;
output.s=s;
output.t=t;
output.f=f;
output.e_f_shaped=e_f_shaped;
output.e_t_shaped=e_t_shaped;
output.params=params;
output.times=times;
output.timePulse=timePulse;
output.nrg=energies;

% *******************************************

% VI: Display Results

% A: Plot State Populations vs. Time
% NOTE: t is array of timesteps from numerical integration
%       times is array of timesteps from pulse shaping

if input.plot==1
    dfigure;
    % Preparing axis limits:
    TI1 = find(times>-T_LIMIT,1);
    TI2 = find(times>T_LIMIT-step,1);

    subplot(2,2,1)
%     plot(t,abs(s).^2);
    plot(times(TI1:TI2),abs(timePulse(TI1:TI2)).^2,'red'); %Plotting Pulse Intensity Profile in Red
    ylabel('Pulse envelope');
    xlabel('Time from Pulse Peak [ps]');
    hold off;

    subplot(2,2,2)
    plot(t,abs(s(:,2)).^2,'green');
    hold on;
    plot(t,abs(s(:,ind2)).^2,'cyan');
    s_ion=s_ion1+s_ion2;
    plot(t,sum(abs(s_ion').^2),'black'); % Plotting Total Ionized Fraction in Black
    legend('|1>','|3>','Ionized');
    titlestr = 'State Populations vs. Time';
%     title(titlestr);
    xlabel('Time from Pulse Peak [ps]');
    ylabel('Population');
    hold off;

    subplot(2,2,4)
    plot(times(TI1:TI2),(abs(timePulse(TI1:TI2)).^2*params.Io/params.I_REF*params.w1+params.w01)/2/pi,'g--');
    hold on;
    plot(times(TI1:TI2),(abs(timePulse(TI1:TI2)).^2*params.Io/params.I_REF*params.w3+params.w03)/2/pi,'c--');
    plot(times(TI1:TI2),(abs(timePulse(TI1:TI2)).^2*params.Io/params.I_REF*params.wC+params.w0C1-params.EMAX/params.HBAR/2*1E-12)/2/pi,'g.');
    plot(times(TI1:TI2),(abs(timePulse(TI1:TI2)).^2*params.Io/params.I_REF*params.wC+params.w0C3-params.EMAX/params.HBAR/2*1E-12)/2/pi,'c.');
    plot(times(TI1:TI2),params.n03*params.wField/2/pi*ones(size(times(TI1:TI2))),'r--');
    plot(times(TI1:TI2),params.n0C3*params.wField/2/pi*ones(size(times(TI1:TI2))),'r--');
    xlabel('Time from Pulse Peak [ps]');
    ylabel('Frequency [THz]')
    titlestr= 'Energy levels and multiphoton transition frequencies';
    title(titlestr);
    legend('|1>','|3>','|C1>','|C3>',[num2str(params.n03) '\cdot h\nu_{field}'],[num2str(params.n0C3) '\cdot h\nu_{field}'])
%     plot(t,abs(s_ion).^2);
%     titlestr = 'Population of Ladder Coefficients vs. Time';


%     ylabel('Population');

    % B: Plot Photoelectron Spectrum
    subplot(2,2,3)
    plot(energies,PHOTOELECTRON_SPECTRUM,'k');
    hold on; plot(energies,output.rsum(:,1),'green')
    hold on; plot(energies,output.rsum(:,2),'cyan')
    hold off
    titlestr = 'Photoelectron Spectrum';
%     title(titlestr);
    xlabel('Energy [eV]');
    ylabel('Photoelectron Yield (Arb.)');
end
    % toc
end
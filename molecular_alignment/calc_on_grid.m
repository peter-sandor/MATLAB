%% Calculation for linear molecules
clear simresult_array
to_scan1=[1:3:31]; % rotational temperatures in [K]
to_scan2=[10:5:50]; % laser intensities in [W/cm2]
Nscan1=length(to_scan1);
Nscan2=length(to_scan2);

ind3=1;
for ind2=1:Nscan2
    for ind1=1:Nscan1
        param_pairs(ind3,:)=[to_scan1(ind1) to_scan2(ind2)];
        ind3=ind3+1;
    end
end
Nscan=Nscan1*Nscan2;
% temp=reshape(param_pairs,[11 11 2]);
%%
for ind1=1:Nscan
    in.Trot=param_pairs(ind1,1); % rotational temperature in [K]
    in.laser_int1=0; % THz peak intensity in [TW/cm2]
    in.laser_int2=param_pairs(ind1,2); % laser peak intensity in [TW/cm2]
    in.t0=0; % relative delay of orientation pulse with respect to alignment pulse
    in.laser_fwhm1=0.05; % THz pulse duration intensity FWHM in [ps]
    in.laser_fwhm2=0.075; % laser pulse duration intensity FWHM in [ps]
    in.polar_antr=31; % polarizability anisotropy [atomic units]
    in.dipole=0.33; % molecular dipole moment [atomic units]
    in.rot_const=0.2026; % rotational constant 'B' in [1/cm]
    in.centrif=1*3.46e-8; % centrifugal distortion 'D' [1/cm]
    in.abund_evenJ=1; % abundance of even J states
    in.abund_oddJ=1; % abundance of odd J states
    in.maxJ=50; % max J level considered
    in.maxdelay=92; % max timedelay in [ps]
    in.timestep=0.1;
    in.Ntheta=200;
    in.rand=0; % random phase for different (J,M) states? (1: yes, 0: no)
    in.solvetype=1; % 0: populate all initial states and solve ('single molecule' calculation); 1: populate a single initial state, solve, then repeat ('infinite number of molecules' approximation)
    in.calc_cos2=1; % calculate cos^2(theta)? (1=yes,0=no)
    in.calc_prob=1; % calculate angular probability distribution? (1=yes,0=no)

    simresult_array(ind1) = solve_align_linear_v2(in);
    
    if mod(ind1,20)==0
        save partial_save.mat simresult_array -v7.3;
    end
    disp([num2str(ind1) '/' num2str(Nscan)]);
end
% input_params=reshape(input_params,[Nscan1 Nscan2]);
simresult_array=reshape(simresult_array,[Nscan1 Nscan2]);
save scan_I_and_T_75fs_OCS_mycode_2.mat simresult_array -v7.3;
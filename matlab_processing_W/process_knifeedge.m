basedir = 'U:\Measurement_Data\spotsize_measurement\2020_06_12_SPECS_spotsize_measurement\scan_A1\'; files1 = dir([basedir '*A2_2820*']);
% basedir = 'U:\Measurement_Data\spotsize_measurement\2020_06_12_SPECS_spotsize_measurement\scan_A2\'; files1 = dir([basedir '*center_2810*']);
N = length(files1);
%%
thrs = 0.03;
ind_fit1=[1 62];

for ind1 = 1:N
    if exist([basedir files1(ind1).name '\data_dot.txt'])==0
        comma2dot([basedir files1(ind1).name '\data.txt']);
    end
    ind_pos = strfind(files1(ind1).name,'A3_')+3;
    N_name = length(files1(ind1).name);
    if N_name-ind_pos>4
        ind_pos2 = strfind(files1(ind1).name,'_center')-1;
    else
        ind_pos2 = N_name;
    end
    str_A3pos = files1(ind1).name(ind_pos:ind_pos2);
    ind_m = strfind(str_A3pos,'m');
    if ~isempty(ind_m)
        pos_A3(ind1) = str2num(['-' str_A3pos(2:end)]);
    else
        pos_A3(ind1) = str2num(str_A3pos);
    end
    data_A1_raw{ind1} = load([basedir files1(ind1).name '\data_dot.txt']);
    N_trace(ind1) = size(data_A1_raw{ind1},1);
    data_A1{ind1}(:,1) = data_A1_raw{ind1}(:,1);
    data_A1{ind1}(:,2) = (data_A1_raw{ind1}(:,2) - min(data_A1_raw{ind1}(:,2)))/(max(data_A1_raw{ind1}(:,2)) - min(data_A1_raw{ind1}(:,2)));
    lineout_A1{ind1}(:,1) = diff(data_A1{ind1}(:,1))/2 + data_A1{ind1}(1:N_trace-1,1);
    lineout_A1{ind1}(:,2) = diff(data_A1{ind1}(:,2))./diff(data_A1{ind1}(:,1))/max(diff(data_A1{ind1}(:,2))./diff(data_A1{ind1}(:,1)));
    lineout_A1{ind1}(:,2) = abs(lineout_A1{ind1}(:,2)/max(lineout_A1{ind1}(:,2)));
    param_A_guess = 1;
    
    if data_A1{ind1}(1,2)<0.5*data_A1{ind1}(end,2)
        [param_C_guess,~,param_B_guess] = peak_props(lineout_A1{ind1});
        param_D_guess = 0;
        FIT0(ind1) = ezfit(data_A1{ind1}(ind_fit1(1):ind_fit1(2),1),data_A1{ind1}(ind_fit1(1):ind_fit1(2),2),['1/2*erf((x-C)/B)+0.5; A = ' num2str(param_A_guess) '; B = ' num2str(param_B_guess) '; C = ' num2str(param_C_guess) '; D = ' num2str(param_D_guess) ';']);
        widths(ind1) = FIT0(ind1).m(1);
    elseif data_A1{ind1}(end,2)<0.5*data_A1{ind1}(1,2)
        [param_C_guess,~,param_B_guess] = peak_props(lineout_A1{ind1});
        param_D_guess = 0;
        FIT0(ind1) = ezfit(data_A1{ind1}(ind_fit1(1):ind_fit1(2),1),data_A1{ind1}(ind_fit1(1):ind_fit1(2),2),['1/2*erfc((x-C)/B)+0.5; A = ' num2str(param_A_guess) '; B = ' num2str(param_B_guess) '; C = ' num2str(param_C_guess) '; D = ' num2str(param_D_guess) ';']);
        widths(ind1) = FIT0(ind1).m(1);
    else
        param_C1_guess0 = mean(lineout_A1{ind1}(logical([(lineout_A1{ind1}(:,2)>0.8).*(lineout_A1{ind1}(:,1)<(lineout_A1{ind1}(1,1)+lineout_A1{ind1}(end,1))/2)]),1));
        param_C2_guess0 = mean(lineout_A1{ind1}(logical([(lineout_A1{ind1}(:,2)>0.8).*(lineout_A1{ind1}(:,1)>(lineout_A1{ind1}(1,1)+lineout_A1{ind1}(end,1))/2)]),1));
        ind_center = value2index(lineout_A1{ind1}(:,1),(param_C1_guess0+param_C2_guess0)/2);
        [param_C1_guess,~,param_B1_guess] = peak_props(lineout_A1{ind1}(1:ind_center,:));
        [param_C2_guess,~,param_B2_guess] = peak_props(lineout_A1{ind1}(ind_center:end,:));
        param_D_guess = -0.5;
        FIT0(ind1) = ezfit(data_A1{ind1}(ind_fit1(1):ind_fit1(2),1),data_A1{ind1}(ind_fit1(1):ind_fit1(2),2),...
        ['A*((1-1/2*erf((x-C1)/B1))+(1-1/2*erfc((x-C2)/B2)))+D; A = ' num2str(param_A_guess) '; B1 = ' num2str(param_B1_guess) '; B2 = ' num2str(param_B2_guess) '; C1 = ' num2str(param_C1_guess) '; C2 = ' num2str(param_C2_guess) '; D = ' num2str(param_D_guess) ';']);
        widths(ind1) = mean(FIT0(ind1).m(2:3));
    end
    
    if 0
        figure;
        subplot(211);hold on;
        plot(data_A1{ind1}(:,1),data_A1{ind1}(:,2),'ko')
        xlim([min(data_A1{ind1}(:,1)),max(data_A1{ind1}(:,1))]);
        showfit(FIT0(ind1))
        subplot(212)
        plot(lineout_A1{ind1}(:,1),lineout_A1{ind1}(:,2),'r')
        xlim([min(data_A1{ind1}(:,1)),max(data_A1{ind1}(:,1))]);
    end
    
    if 0
        y = FIT0(ind1).m(1)/2*erf((x-FIT0(ind1).m(3))/FIT0(ind1).m(2))+FIT0(ind1).m(4);
        figure;plot(x,y);
        lineout_smooth{ind1} = smooth1D(lineout_A1{ind1},1.5);
        x0{ind1} = cumsum(diff(data_A1{ind1}(:,1)));
        x_A1{ind1} = x0{ind1} - x0{ind1}(vec2ind(lineout_smooth{ind1}==max(lineout_smooth{ind1})));
        lineout_calc{ind1} = lineout_smooth{ind1}/max(lineout_smooth{ind1});
        lineout_calc{ind1}(lineout_calc{ind1}<=thrs) = 0;
    end
    
    if 0
        hfig2 = figure;
        subplot(211);plot(data_A1{ind1}(:,2),'k');xlim([1 length(data_A1{ind1}(:,2))]);
        subplot(212);plot(lineout_A1{ind1},'r');xlim([1 length(data_A1{ind1}(:,2))]);
        temp = ginput(2);
        close(hfig2)
        ind_fit(ind1,:) = round(temp(:,1).');
        FIT_A1(ind1) = ezfit(x_A1{ind1}(ind_fit(1,1):ind_fit(1,2)),lineout_calc{ind1}(ind_fit(1,1):ind_fit(1,2)),'exp(-2*(x-A)^2/B^2);A = 0; B = 5;');

        [mean_x(ind1),sigma(ind1)] = peak_props([map2colvec(x_A1{ind1}), map2colvec(lineout_calc{ind1})]);
        FIT_A1(ind1) = ezfit(x_A1{ind1},lineout_calc{ind1},['exp(-2*(x-A)^2/B^2);A = 0; B = ' num2str(sigma(ind1)) ';']);
        figure;
        subplot(211)
        plot(data_A1{ind1}(:,1),data_A1{ind1}(:,2),'k')
        title(num2str(pos_A3(ind1)))
        subplot(212);
        hndl_plot1(ind1) = plot(x_A1{ind1},lineout_calc{ind1},'o--');
        hold on;plot(x_A1{ind1},lineout_A1{ind1}/max(lineout_A1{ind1}),'o--');
        showfit(FIT_A1(ind1));
        widths(ind1) = FIT_A1(ind1).m(2);
    end
end
%%
% legend(hndl_plot1,num2str(map2colvec(pos_A3)));
[pos_sorted,ind_sort] = sort(pos_A3,'ascend');
widths_sorted = widths(ind_sort);
% widths_sorted2 = widths2(ind_sort);
% sigma_sorted = sigma(ind_sort);

figure;plot(pos_sorted,widths_sorted,'ko-')
% hold on;plot(pos_sorted,2*sigma_sorted,'ro-')
xlabel('A3 position [micron]')
ylabel('transition widths [micron]')
legend('from error function fitting')
% legend('from Gaussian fitting','from calculating standard deviation')
setfigP;

if 0
    figure;plot(pos_sorted,widths_sorted2,'ko-')
    hold on;plot(pos_sorted,2*sigma_sorted,'ro-')
    xlabel('A3 position [micron]')
    ylabel('transition widths [micron]')
    legend('from error function fitting')
    legend('from Gaussian fitting','from calculating standard deviation')
    setfigP;
end
%%
ind_fit2 = [1 N];
z = pos_sorted;
w = widths_sorted;
nonlinearity = 2;
lambda = 0.8; % [um]
% FIT_gaussbeam = ezfit(z(ind_fit2(1):ind_fit2(2)),w(ind_fit2(1):ind_fit2(2)),['y=w0*sqrt(1+((x-x0)/(pi*w0^2/' num2str(lambda) '))^2);x0 = 1100; w0=' num2str(min(w))]);
FIT_gaussbeam = ezfit(z,w,['y=M*w0*sqrt(1+((x-x0)/(pi*w0^2/' num2str(lambda) '))^2);x0 = ' num2str(z(w == min(w))) '; w0=' num2str(min(w)) '; M = 1.3;']);
M_fitted = FIT_gaussbeam.m(1);
w0_fitted = FIT_gaussbeam.m(2);
zC_fitted = FIT_gaussbeam.m(3);
w_fitted = M_fitted*w0_fitted*sqrt(1+((z-zC_fitted)/(pi*w0_fitted^2/lambda)).^2);
figure;plot(z,w,'ko-')
% hold on;plot(z,w_fitted,'r-')
xlabel('z [micron]')
ylabel('transition widths [micron]')
setfigP;
showfit(FIT_gaussbeam)
function out = fit2dScan(in)
record=[];
out = in;
% figure;
Np = length(in.N_params);
if in.flag_grid
    for ind1 = 1:Np
        if in.params(ind1,1) == in.params(ind1,2)
            out.vec_param{ind1} = in.params(ind1,1);
        else
            out.vec_param{ind1}=in.params(ind1,1):(in.params(ind1,2)-in.params(ind1,1))/(in.N_params(ind1)-1):in.params(ind1,2);
        end
    end
    out.vec_comb = combine(out.vec_param);
    Ncomb = size(out.vec_comb,1);
    flag1=0;
    temp = Ncomb;
    for ind2 = 1:Ncomb
        out.chi2mat(ind2) = calc_pulse([out.vec_comb(ind2,:) 0]);
        if ind2>1
            if all(out.chi2mat(ind2)<out.chi2mat(1:ind2-1))
                temp = ind2;
            end
        end
        disp([num2str(ind2) '/' num2str(Ncomb)]);
    end
    start_params = out.vec_comb(temp,:);
else
    if ~isempty(in.start_params)
        start_params = in.start_params;
    else
        for ind1 = 1:Np
            start_params(ind1) = (in.params(ind1,2)-in.params(ind1,1))*rand - in.params(ind1,1);
        end
    end
end

str_display='off';
% str_display='iter';
max_fun_eval = 200;
tol_in_x = 1e-2;
tol_in_fun = 1e-2;
if in.flag_plot
    options = optimset('Display',str_display,'PlotFcns',@optimplotfval,'MaxFunEvals',max_fun_eval,'TolX',tol_in_x,'TolFun',tol_in_fun);
else
    options = optimset('Display',str_display,'MaxFunEvals',max_fun_eval,'TolX',tol_in_x,'TolFun',tol_in_fun);
end
flag1=in.flag_progress;
ind_iter=1;
fprintf('%s\n','Performing search')
record=[];
if in.start_param_SPM ~= 0
    [a_out,G_last] = fminsearch(@calc_pulse,[start_params in.start_param_SPM],options);
else
    [a_out,G_last] = fminsearch(@calc_pulse,start_params,options);
end
record_sorted=sortrows(record,size(record,2));
fprintf('%s\n','done.')
disp(in.mask);
disp(record_sorted(1,:));
mask_retrieved = in.mask;
for ind1 = 1:length(in.N_params)
    mask_retrieved = strrep(mask_retrieved,['a(' num2str(ind1) ')'],num2str(record_sorted(1,ind1)));
end

out.record=record_sorted;
out.mask_retrieved = mask_retrieved;
out.E_retrieved = apply_spectral_mask(in.E_om(:,1),in.E_om(:,2),mask_retrieved);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function G_error = calc_pulse(a)
%     mask = ['exp(1i*(' num2str(a(1)) '*1e2*(omega-omega0).^2+' num2str(a(2)) '*1e3*(omega-omega0).^3+'  num2str(a(3)) '*1e4*(omega-omega0).^4))'];
    mask = in.mask;
    for ind1 = 1:Np
        mask = strrep(mask,['a(' num2str(ind1) ')'],num2str(a(ind1)));
    end
    spomega = apply_spectral_mask(in.E_om(:,1),in.E_om(:,2),mask);
    calced_traces = calc_dscan(in.E_om(:,1),spomega,in.material,in.mat_thickness);
    for ind1 = 1:length(in.mat_thickness)
        trace_ip(ind1,:) = interp1(calced_traces.axis_lambda,calced_traces.trace_lambda(ind1,:),in.wl_exp,'spline');
    end
    G_error = calc_G_error(in.trace_exp,trace_ip);
    record=cat(1,record,[map2rowvec(a) G_error]);
    if flag1==1
        fprintf('%s','.')
        if mod(ind_iter,100)==0
            fprintf('%s\n','');
        end
        ind_iter=ind_iter+1;
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
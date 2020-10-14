function output = beam_fitter(input,systems)

minsearch_options.Display='notify';
minsearch_options.MaxFunEvals=1000;
minsearch_options.MaxIter=1000;
minsearch_options.TolFun=1e-2;
minsearch_options.TolX=1e-1;
% minsearch_options.OutputFcn=[];
% minsearch_options.PlotFcns='@optimplotfval';
output.search_path=[];
output.bs_all=[];
[a_out,S_min,exitflag,out0]=fminsearch(@fit_beam,[input.params.w input.params.R input.params.M],minsearch_options);
output.w=a_out(1);
output.R=a_out(2);
output.S_min=S_min;

    function S = fit_beam(a)
    Nsys=length(systems);
    S=0;
    params.w=a(1);
    params.R=a(2);
    params.lambda=input.params.lambda;
    M=a(3);
    for ind1=1:Nsys
        result(ind1)=beam_tracer(params,systems{ind1});
        output.beamsizes(ind1)=M*result(ind1).w;
        S = S + (output.beamsizes(ind1) - input.wmeas(ind1))^2;
    end
    output.search_path=cat(1,output.search_path,[a(1) a(2) a(3) S]);
    output.bs_all = cat(1,output.bs_all,map2rowvec(output.beamsizes));
    end
end
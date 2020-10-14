% function scan_param_peter
clear;
S.pistep=150;
S.GDD=0; % chirp parameter in [fs^2]
S.TOD=0; % chirp parameter in [fs^3]
S.FWHM=0.03; % pulse duration FWHM in [ps]
S.nu_field=390; % central laser frequency in [THz]
S.plot=0; % set this to 1 to get the code display a 4-panel plot for each calculation, or 0 for no plots
scanned=(5:5:50)*1e15; % input parameters to be scanned here.
keep_yield=0;
yield_target=1.0e-3;
yield_tolerance=1e-4;
N=length(scanned);
for ind1=1:N
    ints1=2e16;
    S.intensity=scanned(ind1); % change the name of the variable after the dot to the name of the parameter to be scanned (e.g.: S.intensity, S.FWHM)
%     S.intensity=ints1; % laser intensity in [W/m^2]
    result=SFI_2plus1(S);
    yield1=sum(result.PES(10:end-10));
    if keep_yield
        ints2=2*ints1;
        S.intensity=ints2; % laser intensity in [W/m^2]
        result=SFI_2plus1(S);
        yield2=sum(result.PES(10:end-10));
        while yield2>=(yield_target+yield_tolerance) || yield2<=(yield_target-yield_tolerance)
            ints_new=abs((yield_target-yield1)/(yield2-yield1)*(ints2-ints1)+ints1);
            ints1=ints2;
            ints2=ints_new;
            S.intensity=ints2;
            result=SFI_2plus1(S);
            yield1=yield2;
            yield2=sum(result.PES(10:end-10));
            fprintf('%e, %e\n',ints2,yield2);
        end
    end
    intensity(ind1)=S.intensity;
    PES(:,ind1)=map2colvec(result.PES);
    disp([num2str(ind1) '/' num2str(N)]);
end
%%
% PES_norm=PES./permute(extend(sum(PES,1),size(PES,1)),[2 1]);
% PES_norm2=PES./permute(extend(max(PES,[],1),size(PES,1)),[2 1]);
nrg=map2colvec(result.nrg);
freq=map2colvec(result.f);
to_plot=PES;%[result.rsum(:,1) result.rsum(:,2) map2colvec(result.PES)];
dfigure;

% imagescP(scanned,nrg,to_plot)
% xlabel('\pi-step frequency (relative to center) [THz]')
% ylabel('K [eV]')
% colorbar
% caxis([0 max(max(to_plot(10:end-10,:)))]);

plot(nrg,to_plot)
legend(num2str(map2colvec(scanned)))
xlabel('K [eV]')
setlines;
setfigP;

title(['\Omega_{3C}=' num2str(result.params.rabi3C3/2/pi) ' THz, \Omega_{01}=' num2str(result.params.rabi01/2/pi) ' THz, \chi_{13}=' num2str(result.params.coupling13/2/pi) ' THz'])

% function y=wrapper(x)
%     S.intensity=x;
%     result=SFI_Inverted_A_for_paper(S);
%     yield=sum(result.PES(10:end-10));
%     y=(yield-yield_target)^2;
% end

% end
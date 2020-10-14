function out = fit2IAC(in)

options = optimset('Display','iter','MaxFunEvals',100,'TolX',1e-4,'TolFun',1e-2);
record=[];
% figure;
GDDstart=0.2;
TODstart=0.07;
[a_out,chi2_last] = fminsearch(@calc_pulse,[GDDstart TODstart]);
record0=record;
GDDvec=-1.5*GDDstart:3*GDDstart/9:1.5*GDDstart;
TODvec=-1.5*TODstart:3*TODstart/9:1.5*TODstart;
for ind3=1:length(GDDvec)
    for ind2=1:length(TODvec)
        chi2mat(ind3,ind2) = calc_pulse([GDDvec(ind3) TODvec(ind2)]);
    end
    disp([num2str(ind3) '/' num2str(length(GDDvec))]);
end
figure;imagescP(TODvec,GDDvec,chi2mat)
hold on;plot(record0(:,2),record0(:,1),'r*')

out.record=record;

    function chi2 = calc_pulse(a)
    mask = ['exp(1i*(' num2str(a(1)) '*1e2*(omega-omega0).^2+' num2str(a(2)) '*1e3*(omega-omega0).^3))'];
    [t,Et1,omega,spomega]=simulate_pulse_shaper([in.sp_omega(:,1) in.sp_omega(:,2)],mask);
    tau=in.t_meas;
    Ntau=length(tau);
    [t0 tau0 FWHM]=calc_stats([map2colvec(t) map2colvec(abs(Et1).^2)]);
    clear Scorr;
    for ind1=1:Ntau
        index_start=min(vec2ind(t>=max(t(1),t(1)+tau(ind1))));
        index_end=max(vec2ind(t<=min(t(end),t(end)+tau(ind1))));
        Et2=interp1(t,Et1,t(index_start:index_end)-tau(ind1),'linear');
        Scorr(ind1)=sum(abs((Et1(index_start:index_end)+Et2).^2).^2);
        Scorr2a(ind1)=sum(abs(Et1(index_start:index_end)).^4)+sum(abs(Et2).^4); % constant term
%         Scorr2b(ind1)=4*sum(abs(Et1(index_start:index_end)).^2.*abs(Et2).^2); % intensity autocorrelation
%         Scorr2c(ind1)=4*sum((abs(Et1(index_start:index_end)).^2+abs(Et2).^2).*real(Et1(index_start:index_end).*conj(Et2))); % interferometric term at the fundamental frequency
%         Scorr2d(ind1)=2*sum(real(Et1(index_start:index_end).^2.*conj(Et2).^2)); % interferometric term at the second harmonic
    end
%     Scorr2=Scorr2a+Scorr2b+Scorr2c+Scorr2d;
    index_norm=min(vec2ind(tau>=(tau(1)+2*tau0)));
%     normfactor=Scorr(index_norm);
%     normfactor=max(Scorr)/8;
    normfactor=Scorr2a(round(Ntau/2));
    Scorr=map2colvec(Scorr)/normfactor;
%     Scorr2=Scorr2/normfactor;
    chi2=sum((in.Scorr_meas(index_norm:Ntau-index_norm)-Scorr(index_norm:Ntau-index_norm)).^2);
    record=cat(1,record,[map2rowvec(a) chi2]);
%     plot(tau(index_norm:Ntau-index_norm),in.Scorr_meas(index_norm:Ntau-index_norm),'k')
%     hold on;plot(tau(index_norm:Ntau-index_norm),Scorr(index_norm:Ntau-index_norm),'r')
%     hold off
    end
end
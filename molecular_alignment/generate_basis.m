function out = generate_basis(in)
% fields for input 'in':
%   simresult.data = 2D and 3D angular distributions
%   simresult.delay = timedelay in [ps]
%   simresult.theta = theta values normalized to pi
%   basisnames = cell array with the basis names as strings (e.g. sin(K*theta)^M)
%   K = array with actual values of constant 'K'
%   M = array with actual values of constant 'M'
out=in;

if isfield(in.simresult,'ang_prob')
    calc_3D=squeeze(in.simresult.ang_prob);
    theta=in.simresult.theta;
    delay_calc0 = in.simresult.delay*24.2e-6;
else
    calc_3D=squeeze(in.simresult.data(:,:,1)); % use the '2D' angular distribution (P = |sum_{lm}(a_{lm}*Y_l^m)|^2, without sin(theta) factor)
    theta=pi*in.simresult.theta;
    delay_calc0 = in.simresult.delay;
end
Ndelay=length(delay_calc0);
calc_stepsize=mean(diff(delay_calc0));
Ntheta=length(theta);
theta_step=theta(2)-theta(1);
theta_ext=permute(extend(theta,Ndelay),[2 1]);
calc_3D=calc_3D./extend(trapz(theta,calc_3D.*sin(theta_ext),2),size(calc_3D,2));

basisnames_subs=in.basisnames;
Nbasis=length(in.basisnames);
for indB=1:Nbasis
    basisnames_subs{indB}=StringReplace(basisnames_subs{indB},'K',num2str(in.K(indB)));
    basisnames_subs{indB}=StringReplace(basisnames_subs{indB},'M',num2str(in.M(indB)));
    basisnames_subs{indB}=StringReplace(basisnames_subs{indB},'\theta','theta_ext');
    basisnames_subs{indB}=StringReplace(basisnames_subs{indB},'^','.^');
end

ionRate=zeros([Ndelay Ntheta Nbasis]);
for indB=1:Nbasis
    ionRate(:,:,indB) = eval(basisnames_subs{indB});
    basis(:,indB)=trapz(theta,ionRate(:,:,indB).*calc_3D.*sin(theta_ext),2);
end
out.basis=basis;
end
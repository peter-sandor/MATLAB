function plotreconst2(varargin)

temp = dir('*Ek*');
Ek=load(temp(1).name);
temp = dir('*Speck*');
Speck=load(temp(1).name);
% temp = dir('*frog*');
% fid=fopen(temp(1).name);
% settings = textscan(fid, '%s');
% fclose(fid);
% for ind1=1:length(settings{1})
%     if strcmp('error',settings{1}{ind1})
%     str_frogerror=cell2mat([settings{1}(ind1-1) ' ' settings{1}(ind1) ' ' settings{1}(ind1+1) ' ' num2str(roundP(str2num(settings{1}{ind1+2})*100,2)) '%']);
%     end
% end

reconst_omega=SpectrumConvert(Speck(:,1:3));
[a b tau temp crop]=calc_stats([map2colvec(1:size(Ek,1)) Ek(:,2)]);
cropt_start=crop(1);
cropt_end=crop(2);

[a b tau temp crop]=calc_stats([map2colvec(1:size(Speck,1)) Speck(:,2)]);
crops_start=crop(1);
crops_end=crop(2);

[a b tau temp]=calc_stats(Ek(cropt_start:cropt_end,1:2));
textcoord=temp;

tphase=unwrap(angle(Ek(:,4)+i*Ek(:,5)));
sphase=unwrap(angle(Speck(:,4)+i*Speck(:,5)));

if nargin==1
    spectrum_in=spproc2(varargin{1}); %load an experimentally measured spectrum
    [measured_omega,measured_lambda]=SpectrumConvert(spectrum_in); % omega,E(omega) and lambda,E(lambda)
    sphase_ip=interp1(reconst_omega(:,1),reconst_omega(:,3),measured_omega(:,1),'linear');
    if min(measured_lambda(:,1))<min(Speck(crops_start:crops_end,1))
        crops_end=max((Speck(:,1)>min(measured_lambda(:,1))).*map2colvec(1:size(Speck(:,1))));
    end
    if max(measured_lambda(:,1))>max(Speck(crops_start:crops_end,1))
        indvec1=unique((Speck(:,1)>max(measured_lambda(:,1))).*map2colvec(1:size(Speck(:,1))));
        indvec1(indvec1==0)=[];
        crops_start=max(indvec1);
    end
    measured_spectrum_ip=interp1(measured_lambda(:,1),measured_lambda(:,2).^2,Speck(crops_start:crops_end,1),'linear');
    td_field_from_spectrum=simulate_pulse_shaper([measured_omega(:,1) measured_omega(:,2).*exp(i*sphase_ip)],'1');
    td_ip=interp1(-td_field_from_spectrum(:,1),td_field_from_spectrum(:,2),Ek(cropt_start:cropt_end,1),'linear');
else
    sphase_ip=[];
    td_ip=[];
    measured_spectrum_ip=[];
end

hndl1=figure;
subplot(211)
text_td{1}='Delay [fs]';
text_td{2}='Intensity (normalized)';
text_td{3}='Phase [\pi]';
text_td{4}=[];
% text_td{4}=str_frogerror;
plotyyP2(Ek(cropt_start:cropt_end,1),[Ek(cropt_start:cropt_end,2) abs(td_ip).^2],tphase(cropt_start:cropt_end)/pi,text_td)
% plotyyP2(Ek(cropt_start:cropt_end,1),[Ek(cropt_start:cropt_end,2) abs(td_ip).^2],tphase(cropt_start:cropt_end)/pi)
text(textcoord(1),textcoord(2),[num2str(roundP(tau,1)),' fs'],'Fontsize', 14)

subplot(212)
text_sd{1}='Wavelength [nm]';
text_sd{2}='Spectral density [arb. units]';
text_sd{3}='Phase [\pi]';
text_sd{4}=[];
plotyyP2(Speck(crops_start:crops_end,1),[Speck(crops_start:crops_end,2) measured_spectrum_ip],sphase(crops_start:crops_end)/pi,text_sd)

saveas(hndl1,'reconstructed.fig')
close(hndl1)
end
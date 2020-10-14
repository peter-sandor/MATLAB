function plotreconst(varargin)

if nargin==0
    Ek=load('Ek.dat');
    Speck=load('Speck.dat');
elseif nargin==1
    Ek=load(strcat(varargin{1},'Ek.dat'));
    Speck=load(strcat(varargin{1},'Speck.dat'));
end

hndl2=figure;

plot(Ek(:,2),'k')
title('Select region of interest in time domain')
temp=ginput(2);
cropt_start=uint16(round(min(temp(1,1),temp(2,1))));
cropt_end=uint16(round(max(temp(1,1),temp(2,1))));

plot(Speck(:,2),'k')
title('Select region of interest in spectral domain','color','k')
temp=ginput(2);
crops_start=uint16(round(min(temp(1,1),temp(2,1))));
crops_end=uint16(round(max(temp(1,1),temp(2,1))));

plot(Ek(cropt_start:cropt_end,1), Ek(cropt_start:cropt_end,2),'k')
hold on;plot(Ek(cropt_start:cropt_end,1),0.5*ones([cropt_start-cropt_end+1,1]),'r')
title('Select points for FWHM')
temp=ginput(2);
textcoord=temp(2,:);
tau=abs(diff(temp(:,1)));
close(hndl2)

tphase=unwrap(angle(Ek(:,4)+i*Ek(:,5)))/pi;
sphase=unwrap(angle(Speck(:,4)+i*Speck(:,5)))/pi;

hndl1=figure;
subplot(121)
[AX, H1, H2] = plotyy (Ek(cropt_start:cropt_end,1), Ek(cropt_start:cropt_end,2), Ek(cropt_start:cropt_end,1), tphase(cropt_start:cropt_end));
set (get(AX(1), 'Ylabel'), 'String', 'Intensity [arb. units]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'blue');
set (get(AX(2), 'Ylabel'), 'String', 'Phase [\pi]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'magenta');
set (AX(1), 'Ycolor', 'blue','Ylim',[min(Ek(cropt_start:cropt_end,2)),max(Ek(cropt_start:cropt_end,2))],'Ytick',[min(Ek(cropt_start:cropt_end,2)):(max(Ek(cropt_start:cropt_end,2))-min(Ek(cropt_start:cropt_end,2)))/4:max(Ek(cropt_start:cropt_end,2))]);
set (AX(2), 'Ycolor', 'magenta','Ylim',[min(tphase(cropt_start:cropt_end)),max(tphase(cropt_start:cropt_end))],'Ytick',[min(tphase(cropt_start:cropt_end)):(max(tphase(cropt_start:cropt_end))-min(tphase(cropt_start:cropt_end)))/4:max(tphase(cropt_start:cropt_end))]);
set(H2, 'Linestyle', '-', 'Linewidth', 2, 'Color', 'magenta');
set(H1, 'Linestyle', '-', 'Linewidth', 2, 'Color', 'blue');
% set(AX,'xlim',[nu_lim1, nu_lim2]);
%xlim ([nu_lim1 nu_lim2]);
xlabel(AX(1),'Delay [fs]', 'fontname', 'verdana', 'fontsize', 14);
text(textcoord(1),textcoord(2),[num2str(tau),' fs;'],'Fontsize', 14)

subplot(122)
[AX, H1, H2] = plotyy (Speck(crops_start:crops_end,1), Speck(crops_start:crops_end,2), Speck(crops_start:crops_end,1), sphase(crops_start:crops_end));
set (get(AX(1), 'Ylabel'), 'String', 'Spectral density [arb. units]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'blue');
set (get(AX(2), 'Ylabel'), 'String', 'Phase [\pi]', 'fontname', 'verdana', 'fontsize', 14, 'color', 'magenta');
set (AX(1), 'Ycolor', 'blue','Ylim',[min(Speck(crops_start:crops_end,2)),max(Speck(crops_start:crops_end,2))],'Ytick',[min(Speck(crops_start:crops_end,2)):(max(Speck(crops_start:crops_end,2))-min(Speck(crops_start:crops_end,2)))/4:max(Speck(crops_start:crops_end,2))]);
set (AX(2), 'Ycolor', 'magenta','Ylim',[min(sphase(crops_start:crops_end)),max(sphase(crops_start:crops_end))],'Ytick',[min(sphase(crops_start:crops_end)):(max(sphase(crops_start:crops_end))-min(sphase(crops_start:crops_end)))/4:max(sphase(crops_start:crops_end))]);
set(H2, 'Linestyle', '-', 'Linewidth', 2, 'Color', 'magenta');
set(H1, 'Linestyle', '-', 'Linewidth', 2, 'Color', 'blue');
% set(AX,'xlim',[nu_lim1, nu_lim2]);
%xlim ([nu_lim1 nu_lim2]);
xlabel(AX(1),'Wavelength [nm]', 'fontname', 'verdana', 'fontsize', 14);
% text(lambda(1)*1e9,(.1*(max(rld)-min(rld))+min(rld))*1e6,[num2str(G),' groove/mm;',' f = ' num2str(f), ' mm;',' spectral order: ',num2str(m)], 'Fontsize', 14)

if nargin==0
    saveas(hndl1,'reconstructed.fig')
elseif nargin==1
    saveas(hndl1,strcat(varargin{1},'reconstructed2.fig'))
end
end
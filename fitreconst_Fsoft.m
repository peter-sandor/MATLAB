function varargout = fitreconst(varargin)

if nargin==0
    Speck=load('Speck.dat');
elseif nargin==1
    Speck=load(varargin{1});
end
sphase=unwrap(angle(Speck(:,4)+i*Speck(:,5)));

[speckmean,speckSigma,speckFWHM,indFWHM,ind_nonzero]=calc_stats(Speck(:,1:2));

hndl2=figure;
subplot(211)
plot(sphase(ind_nonzero(1):ind_nonzero(2)),'k')
subplot(212)
plot(Speck(ind_nonzero(1):ind_nonzero(2),2),'k')
% hold on;plot(sphase,'r')
title('Select region of interest in spectral domain','color','k')
temp=ginput(2);
crops_start=uint16(ind_nonzero(1)+round(min(temp(1,1),temp(2,1))));
crops_end=uint16(ind_nonzero(1)+round(max(temp(1,1),temp(2,1))));
close(hndl2)

omega=2*pi*300./Speck(crops_start:crops_end,1);
spectrum=Speck(crops_start:crops_end,2)./omega.^2;
spectrum=spectrum/max(spectrum);
sphase=sphase(crops_start:crops_end)-sphase(round((crops_end+crops_start)/2));
omegamean=sum(spectrum.*omega)/sum(spectrum);
dataout=[omega spectrum sphase];

answer=inputdlg('What orders do you want to fit? (Separate by commas, e.g.: 1,3,4)');
orders=str2num(answer{1});
fittext='y=A0';
% guesstext=[';x0=' num2str(omegamean) ';A0=' num2str(interp1(omega,sphase,omegamean)) ';'];
guesstext=[';A0=' num2str(interp1(omega,sphase,omegamean)) ';'];
for ind1=1:length(orders)
    fittext=[fittext '+A' num2str(orders(ind1)) '*(x-' num2str(omegamean) ')^' num2str(orders(ind1))];
%     answer=inputdlg(['initial guess for order ' num2str(orders(ind1))]);
%     guess(ind1)=str2num(answer{1});
%     guesstext=[guesstext 'A' num2str(orders(ind1)) '=' answer{1} ';'];
    guess(ind1)=10^orders(ind1);
    guesstext=[guesstext 'A' num2str(orders(ind1)) '=' num2str(guess(ind1)) ';'];
end
FIT=ezfit(omega,sphase,[fittext guesstext]);
% out2(ind1,ind2,ind3)=eval_fit_align(in2);
% [Pval(ind1,ind2,ind3-1),Fval(ind1,ind2,ind3-1)] = F_test(N,ind3-1,ind3,chi2(ind1,ind2,ind3-1),chi2(ind1,ind2,ind3));

figure;
subplot(211)
plot(omega,spectrum,'k');
xlabel('\omega [rad/fs]')
ylabel('Intensity [normalized]')
setfigP;

subplot(212)
plot(omega,sphase,'k');
xlabel('\omega [rad/fs]')
ylabel('Phase [rad]')
setfigP
showfit(FIT);

if nargout>=1
    varargout{1}=FIT;
    if nargout==2
        varargout{2}=dataout;
    end
end
end
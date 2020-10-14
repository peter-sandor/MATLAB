load M:\data5\2013_04_30_R\CS2_set3_proc\data.mat
load M:\data5\2013_04_30_R\CS2_set3_proc\pnrg.mat
nrg=data(1).nrg/1000;
Elims=[1.1 2.5];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=2:5;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A1=350;A2=2;A3=0.1;A4=0;';
%%
load M:\data5\2013_04_30_R\CS2_set1_proc\data.mat
load M:\data5\2013_04_30_R\CS2_set1_proc\pnrg.mat
nrg=data(1).nrg/1000;
Elims=[1.1 2.4];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=1:4;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A1=350;A2=2;A3=0.1;A4=0;';
%%
load M:\data5\2013_10_11_R\scan3_proc\data.mat
load M:\data5\2013_10_11_R\scan3_proc\pnrg.mat
nrg=data(1).nrg/1000;
Elims=[1.2 2.7];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=2:6;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A=400;A2=2;A3=0.1;A4=80;';
%%
load M:\data5\2013_10_11_R\scan2_proc\data.mat
load M:\data5\2013_10_11_R\scan2_proc\pnrg.mat
nrg=data(1).nrg/1000;
% Elims=[1.2 2.7];
Elims=[0 1.2];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=1:11;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A1=350;A2=0.6;A3=0.1;A4=0;';
%%
load M:\data5\2013_10_11_R\scan2_proc\data.mat
load M:\data5\2013_10_11_R\scan2_proc\pnrg.mat
nrg=data(1).nrg/1000;
% Elims=[1.2 2.7];
Elims=[0 1.2];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=1:11;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A1=350;A2=0.6;A3=0.1;A4=0;';
%%
load N:\2014_09_19_R\scan17_proc\data.mat
load N:\2014_09_19_R\scan17_proc\pnrg.mat
nrg=data(1).nrg2/1000;
% Elims=[1.2 2.7];
Elims=[0 1.2];
Einds=[value2index(nrg,Elims(1)) value2index(nrg,Elims(2))];
fit_inds=1:4;
fitstring='y=A1*exp(-(x-A2)^2/A3^2)+A4;A1=50;A2=0.8;A3=0.1;A4=0;';
%%
dfigure;
hold on;
for ind1=1:length(data)
    nrg_temp=nrg(Einds(1):Einds(2));
    temp=data(ind1).lnte2(Einds(1):Einds(2),3);
    nrg_temp(temp<=0)=[];
    nrg_tofit{ind1}=nrg_temp;
    temp(temp<=0)=[];
    lnte_tofit{ind1}=temp;
    weights=[];
%     weights=map2rowvec(1./sqrt(lnte_tofit{ind1}));%/sum(1./sqrt(lnte_tofit{ind1}));
%     weights=1/length(temp)*ones([1 length(temp)]); uniform weighing
    FIT_temp=ezfit(nrg_tofit{ind1},[map2rowvec(lnte_tofit{ind1}); weights],fitstring);
    FIT(ind1)=eval_fit(ezfit(nrg_tofit{ind1},[map2rowvec(lnte_tofit{ind1}); map2rowvec(sqrt(lnte_tofit{ind1}))],fitstring));
    peak_pos(ind1)=FIT(ind1).m(2);
%     peak_pos_sigma(ind1)=FIT(ind1).B(2,2);
%     peak_pos2(ind1)=calc_stats([nrg_tofit{ind1} lnte_tofit{ind1}-FIT(ind1).m(4)]);
    
%     errorbar(FIT(1).x,FIT(1).y,FIT(1).dy,'k')
    plot(FIT(ind1).x,FIT(ind1).y,'k')
%     str_legend{ind1}=num2str(roundP(pnrg(ind1),2));
    showfit(FIT(ind1));
    str_legend{2*ind1-1}=num2str(roundP(pnrg(ind1),2));
    str_legend{2*ind1}='fit';
end
hold off;
xlabel('K [eV]')
xlim([nrg(Einds(1)) nrg(Einds(2))])
legend(str_legend)
setfigP
%% Just fly in full manual
peak_pos=[834 732 651 523 428 428 428 428 428]/1000; 

%%
% FIT2=eval_fit(ezfit(pnrg(fit_inds),[map2rowvec(peak_pos2(fit_inds)); map2rowvec(peak_pos_sigma(fit_inds))],'y=A*x+B;A=50;B=25000;'));
FIT2=ezfit(pnrg(fit_inds),peak_pos(fit_inds),'y=A*x+B;A=50;B=25000;');
figure;
% errorbar(pnrg,peak_pos,peak_pos_sigma,'k');
% hold on;
plot(pnrg(fit_inds),peak_pos(fit_inds),'m+');
xlabel('pulse energy [\muJ]')
ylabel('peak position [eV]')
setfigP

showfit(FIT2,'fitcolor','k','fitlinestyle','--')
% showfit(FIT2(2),'fitcolor','m','fitlinestyle','--')
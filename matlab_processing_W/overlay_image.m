%%
% pos1 = [-2350 -2150 3; -70 110 3; 525 525 1];
% temp=generate_AttoCube_positions(pos1,'pos_test.txt');
%%
pic_SEM0 = monopic(imread('U:\Nanostructures\Z1_minta\2019-11-29_SEM_Z1\A7-Global.png'));
[pic_eScan0 pic_ax] = construct_image_from_scan('scan2\data_dot.txt');

pic_SEM = flipud(pic_SEM0);
ind_H_SEM = [100 3950];
ind_V_SEM = [647 4485];

pic_eScan = pic_eScan0;
ind_H_eScan = [1 49];
ind_V_eScan = [1 49];

figure;imagescP(pic_SEM); colormap gray;
figure;imagescP(pic_eScan); colormap jet;
% um_per_pixel_SEM = 20/757;
%%
pic_SEM0 = monopic(imread('U:\Nanostructures\Z1_minta\2019-11-29_SEM_Z1\A1-Global.png'));
[pic_eScan0 pic_ax] = construct_image_from_scan('scan3\data_dot.txt');

pic_SEM = flipud(pic_SEM0);
ind_H_SEM = [126 3994];
ind_V_SEM = [560 4355];

pic_eScan = pic_eScan0;
ind_H_eScan = [11 62];
ind_V_eScan = [6 51];

figure;imagescP(pic_SEM); colormap gray;
figure;imagescP(pic_eScan); colormap jet;
% um_per_pixel_SEM = 20/757;
%%
pic_SEM0 = monopic(imread('U:\Nanostructures\Z1_minta\2019-11-29_SEM_Z1\A2-Global.png'));
[pic_eScan0 pic_ax] = construct_image_from_scan('scan4\data_dot.txt');

pic_SEM = flipud(pic_SEM0);
ind_H_SEM = [200 4050];
ind_V_SEM = [800 4480];

pic_eScan = pic_eScan0;
ind_H_eScan = [1 62];
ind_V_eScan = [1 59];

figure;imagescP(pic_SEM); colormap gray;
figure;imagescP(pic_eScan); colormap jet;
% um_per_pixel_SEM = 20/757;
%%
pic_SEM0 = monopic(imread('U:\Nanostructures\Z1_minta\2019-11-29_SEM_Z1\A4-Global.png'));
[pic_eScan0 pic_ax] = construct_image_from_scan('scan5\data_dot.txt');

pic_SEM = flipud(pic_SEM0);
ind_H_SEM = [6 4096];
ind_V_SEM = [516 4566];

pic_eScan = pic_eScan0;
ind_H_eScan = [5 75];
ind_V_eScan = [1 71];

figure;imagescP(pic_SEM); colormap gray;
figure;imagescP(pic_eScan); colormap jet;
% um_per_pixel_SEM = 20/757;
%%
figure;
h1 = axes('Position',[0 0 1 1]);p1=imagescP(pic_SEM);
set(h1,'ydir','normal','XTick',[],'YTick',[],'dataaspectratio',[1 1 1]);
colormap(h1,'gray');
% Foreground image
% figure;
h2=axes('Position',[0 0 1 1]);
p2=imagescP(pic_eScan(ind_V_eScan(1):ind_V_eScan(2),ind_H_eScan(1):ind_H_eScan(2)),'Xdata',ind_H_SEM,'Ydata',ind_V_SEM,'AlphaData',0.3);
colormap(h2,'jet');
set(h2,'color','none','visible','off','ydir','normal','XTick',[],'YTick',[],'dataaspectratio',[1 1 1])
linkaxes([h1 h2])


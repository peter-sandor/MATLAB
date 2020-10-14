comma2dot('data.txt');
[pic_eScan0 pic_ax] = construct_image_from_scan('data_dot.txt');
%%
hfig1 = figure;
hax1 = imagescP(pic_ax{1},pic_ax{2},pic_eScan0);
xlabel('A1 position [\mum]')
ylabel('A2 position [\mum]')
title('electron signal from sample ???')
% set(hax1,'gridalpha',0.5,'dataaspectratio',[1 1 1])
grid on;
colorbar;
colormap jet;
setfigP;
saveas(hfig1,'scan_image.fig')
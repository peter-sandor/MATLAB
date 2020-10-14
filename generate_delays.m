delay_vec=map2colvec([-8:-1 -0.75:0.25:1.5 2:2:30 30.5:0.5:35 35.25:0.25:46 46.5:0.5:50 51:2:70 70.5:0.5:76 76.25:0.25:88 88.5:0.5:92])*1000;
delay_vec=flipdim(delay_vec,1);
dlmwrite('delays_for_scan.txt',delay_vec,'\t');
%%
delay_vec2=map2colvec([-500:40:500]);
delay_vec2=flipdim(delay_vec2,1);
dlmwrite('delays_for_scan2.txt',delay_vec2,'\t');
%%
delay_vec2=map2colvec([72000:250:92000]);
delay_vec2=flipdim(delay_vec2,1);
dlmwrite('delays_for_scan3.txt',delay_vec2,'\t');
%%
delay_vec=map2colvec([-8:-1 -0.75:0.25:1.5 2:2:18 18.5 19:0.25:22 22.5 23 23.5:2:29.5 30.5:0.5:35 35.25:0.25:46 46.5:0.5:50 51:2:58 58.5 59.5:0.25:62 62.5 63.5:2:69.5 70.5:0.5:76 76.25:0.25:88 88.5:0.5:92])*1000;
delay_vec=flipdim(delay_vec,1);
dlmwrite('m8ps_to_92ps_nonequiv.txt',delay_vec,'\t');
%%
delay_vec=map2colvec([-8:-1 -0.5 0 0.5 2:2:18 18.5 19:0.25:22 22.5 23 23.5:2:29.5 30.5:0.5:35 35.25:0.25:46 46.5:0.5:50])*1000;
delay_vec=flipdim(delay_vec,1);
dlmwrite('m8ps_to_p50ps_nonequiv.txt',delay_vec,'\t');
%%
delay_vec=map2colvec([-8:-1 -0.5 0 0.5 2:2:18 18.5 19:0.25:22 22.5 23 23.5:2:29.5 30.5:0.5:35 35.25:0.25:46 46.5:0.5:50 51:2:58 58.5 59.5:0.25:62 62.5 63.5:2:69.5 70.5:0.5:76 76.25:0.25:88 88.5:0.5:92])*1000;
delay_vec=flipdim(delay_vec,1);
dlmwrite('m8ps_to_92ps_nonequiv_v2.txt',delay_vec,'\t');
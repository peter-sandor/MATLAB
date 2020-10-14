%% run beam_fitter for seed beam in amplifier 3rd stage
input.params.w=1;
input.params.R=1e4;
input.params.lambda=8e-4; % mm
input.params.M=1.6;
wx = [2.1 4 3.3 0.66 0.77 0.5 0.66];
wy = [1.67 3.27 3 0.83 0.55 0.5 0.48];
foldername='I:\data2018\2018_04_10_PLSB1\files_for_fit\';
%% run beam fitter for Photonics laser beam in amplifier 3rd stage
input.params.w=0.65;
input.params.R=1e4;
input.params.lambda=5.27e-4; % mm
input.params.M=3.29;
wx = [3.555 7.013 12.666 8.598 10.873 14.627 9.550];
wy = [2.518 5.314 9.29 6.21 7.87 11.67 7.54];
foldername='I:\data2018\2018_04_12_PLSB1\files_for_fit\';
%% run beam fitter for Photonics laser beam in amplifier 3rd stage
input.params.w=0.65;
input.params.R=1e4;
input.params.lambda=5.27e-4; % mm
input.params.M=3.29;
wx = [3.55 7.088 13.366 8.397 10.2 17.745 9.997];
wy = [2.735 5.69 10.095 6.551 8.041 11.379 7.403];
foldername='I:\data2018\2018_04_13_PLSB1\files_for_fit\';
%%
input.wmeas = mean([wx; wy],1);
% input.wmeas = wx;
files_to_read{1}=[foldername 'P1.txt'];
files_to_read{2}=[foldername 'P2.txt'];
files_to_read{3}=[foldername 'P3.txt'];
files_to_read{4}=[foldername 'P4.txt'];
files_to_read{5}=[foldername 'P5.txt'];
files_to_read{6}=[foldername 'P6.txt'];
files_to_read{7}=[foldername 'P7.txt'];

for indF=1:length(files_to_read)
    fid=fopen(files_to_read{indF});
    frewind(fid);
%     system(indF)=[];
    ind1=1;
    while ~feof(fid)
        temp=fgetl(fid);
        if temp~=-1
            systems{indF}(ind1).element=temp(1:3);
            systems{indF}(ind1).param=str2num(temp(4:end));
        end
        ind1=ind1+1;
    end
end
output = beam_fitter(input,systems);
output.bs_all(end,:)./input.wmeas
output.search_path(end,:)
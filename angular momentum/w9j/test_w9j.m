%values taken from Varshoalovich; Quantum Theory of Angular Momentum.
%(1988). Table 10.13

clear all; close all; clc;

c = [0 1 2];
f = [0 1 2];

a = {[0 1 1 2 2 3 3 4],...
    [1 1 0 1 1 2 0 1 1 2 2 1 1 2 2 3 1 2 2 3 3 2 2 3 3 4 2 3 3 4 3 3 4],...
    [2 2 1 2 2 3 1 1 2 2 3 3 0 1 1 2 2 3 3 4 0 1 1 2 2 3 3 4 1 1 2 2 3 3 4 1 2 2 3 3 4 2 2 3 3 4]};
d = {[1 1 3 3 5 5 7 7]/2,...
    [1 3 1 1 3 3 1 1 3 3 5 1 3 3 5 5 3 3 5 5 7 3 5 5 7 7 5 5 7 7 5 7 7]/2,...
    [3 5 3 3 5 5 1 3 3 5 5 7 1 1 3 3 5 5 7 7 1 1 3 3 5 5 7 7 1 3 3 5 5 7 7 3 3 5 5 7 7 3 5 5 7 7]/2};
b = {[0 1 1 2 2 3 3 4],...
    [0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4],...
    [0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4]};
e = {[1 1 3 3 5 5 7 7]/2,...
    [1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 5 5 5 5 7 7 7 7 7 7 7]/2,...
    [1 1 1 1 1 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7 7]/2};

kk = 1;
for ii = 1:3
   for jj = 1:length(a{ii})
       W(kk) = w9j(a{ii}(jj),b{ii}(jj),c(ii),d{ii}(jj),e{ii}(jj),f(ii),1/2,1/2,0);
       kk = kk+1;
   end
end

table1 = [a{1}' d{1}'*2 b{1}' e{1}'*2 W(1:8)'];
table2 = [a{2}' d{2}'*2 b{2}' e{2}'*2 W(9:8+33)'];
table3 = [a{3}' d{3}'*2 b{3}' e{3}'*2 W(8+33+1:87)'];

    

%---print values to text file to print and compare with table values---%
[filename,pathname] = uigetfile;
fid = fopen([pathname,filename], 'w');
fprintf(fid,'%i %i/2 %i %i/2 %+12.6f \r\n',table1');
fprintf(fid,' \r\n')
fprintf(fid,'%i %i/2 %i %i/2 %+12.6f \r\n',table2');
fprintf(fid,' \r\n')
fprintf(fid,'%i %i/2 %i %i/2 %+12.6f \r\n',table3');
fclose(fid);
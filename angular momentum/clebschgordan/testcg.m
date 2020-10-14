clear all; close all; clc

[filename,pathname] = uigetfile('*.txt');

fid = fopen([pathname,filename]);
A = textscan(fid,'%*c %s %s %*c %s %s %*c %s %s %f','CommentStyle',{'=','='});
fclose(fid);

j1 = zeros(7392,1);
j2 = zeros(7392,1);
m1 = zeros(7392,1);
m2 = zeros(7392,1);
j = zeros(7392,1);
m = zeros(7392,1);
C = zeros(7392,1);
for i = 1:7392
    j1(i) = eval(A{1}{i});
    j2(i) = eval(A{2}{i});
    m1(i) = eval(A{3}{i});
    m2(i) = eval(A{4}{i});
    j(i) = eval(A{5}{i});
    m(i) = eval(A{6}{i});
    C(i) = A{7}(i);
end

Cnew = zeros(7392,1);
for i = 1:7392
    Cnew(i) = clebschgordan(j1(i),m1(i),j2(i),m2(i),j(i),m(i));
end

hist(Cnew-C,100)
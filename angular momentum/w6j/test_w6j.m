clear all; close all; clc

[filename,pathname] = uigetfile('*.txt');

fid = fopen([pathname,filename]);
A = textscan(fid,'%*c %s %s %s %*c %s %s %s %s','CommentStyle',{'=','='});
fclose(fid);

a = zeros(2264,1);
b = zeros(2264,1);
c = zeros(2264,1);
d = zeros(2264,1);
e = zeros(2264,1);
f = zeros(2264,1);
W = zeros(2264,1);
for i = 1:2264
    a(i) = eval(A{1}{i});
    b(i) = eval(A{2}{i});
    c(i) = eval(A{3}{i});
    d(i) = eval(A{4}{i});
    e(i) = eval(A{5}{i});
    f(i) = eval(A{6}{i});
    W(i) = eval(A{7}{i});
end

Wnew = zeros(2264,1);
for i = 1:2264
    Wnew(i) = w6j(a(i),b(i),c(i),d(i),e(i),f(i));
end

hist(Wnew-W,100)
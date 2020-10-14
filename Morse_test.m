x=0:0.001:1;
x0=0.2;
D=10;
khsi=0.1;
n=0;
[a,b,c]=morse(x,x0,khsi,D,n);
dfigure;hold off;plot(x,a,'k-');
hold on;plot(x,c+permute(extend(b,length(x)),[1 2]))
ylim([-D 0])

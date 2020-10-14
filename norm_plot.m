function norm_plot(varargin)

if nargin==0
    hndl=gca;
else
    hndl=varargin{1};
end
plotdata=getplot(hndl);
Ny=length(plotdata.y);
for ind1=1:Ny
    ymaxes(ind1)=max(plotdata.y{ind1});
end
ymaxx=max(ymaxes);
for ind1=1:Ny
    plotdata.y{ind1}=plotdata.y{ind1}/ymaxx;
end
plotdata.ylim=[0 1];
figure;
putplot(plotdata);
setfigP;
setlines;
end
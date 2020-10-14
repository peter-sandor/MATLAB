function putplot(struct_plot,varargin)

if nargin==1
	plottype='default';
elseif nargin==2
    plottype=varargin;
end
hndl2=gca;
if iscell(struct_plot.x)
    N=length(struct_plot.x);
    hold on
    if strcmp(plottype,'bar')
        to_plot=[];
        for ind1=1:N
            to_plot = [to_plot map2colvec(struct_plot.y{ind1})];
        end
        bar(hndl2,struct_plot.x{ind1},to_plot);
    else
        for ind1=1:N
            plot(hndl2,struct_plot.x{ind1},struct_plot.y{ind1});
        end
    end
    hold off
    if ~isempty(struct_plot.z{1})
        legend(struct_plot.z{1})
    end
else
    M=size(struct_plot.z);
    N=1;
    if (M(1)*M(2)==M(1)) || (M(1)*M(2)==M(2))
        if strcmp(plottype,'bar')
            bar(hndl2,struct_plot.x,struct_plot.y);
        else
            plot(hndl2,struct_plot.x,struct_plot.y);
        end
    elseif M(1)>1 && M(2)>1
        if strcmp(plottype,'default')
            pcolor(struct_plot.x,struct_plot.y,struct_plot.z)
        elseif strcmp(plottype,'imagesc')
            imagesc(struct_plot.x,struct_plot.y,struct_plot.z)
        elseif strcmp(plottype,'imagescP')
            imagescP(struct_plot.x,struct_plot.y,struct_plot.z)
        elseif strcmp(plottype,'surf')
            surf(struct_plot.x,struct_plot.y,struct_plot.z)
            view([0 90]);
            shading interp;
        elseif strcmp(plottype,'pcolor')
            pcolor(struct_plot.x,struct_plot.y,struct_plot.z)
        elseif strcmp(plottype,'contour')
            contour(struct_plot.x,struct_plot.y,struct_plot.z)
        end
        caxis(struct_plot.clim);
    else
        plot(hndl2,struct_plot.x,struct_plot.y)
    end
end

xlim(struct_plot.xlim);
ylim(struct_plot.ylim);
if isfield(struct_plot,'title')
    title(hndl2,struct_plot.title)
end
if isfield(struct_plot,'xlabel')
    xlabel(hndl2,struct_plot.xlabel)
end
if isfield(struct_plot,'ylabel')
    ylabel(hndl2,struct_plot.ylabel)
end
if isfield(struct_plot,'xscale')
    set(hndl2,'xscale',struct_plot.xscale)
end
if isfield(struct_plot,'yscale')
    set(hndl2,'yscale',struct_plot.yscale)
end
if N>1 && isempty(struct_plot.z)==0
    legend(hndl2,struct_plot.z{1})
end
% setlines
% setfigP
end
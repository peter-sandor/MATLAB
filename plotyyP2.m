function plotyyP2(varargin)
% function plotyyP2(x,y1,y2,texts)

if nargin>=3
    x=varargin{1};
    y1=varargin{2};
    y2=varargin{3};
    y=[y1 y2];
    for ind1=1:2
        if min(y(:,ind1))<0;
            factormin(ind1)=1.0;%1.1;
        else factormin(ind1)=1.0;%0.9;
        end
        if max(y(:,ind1))<0;
            factormax(ind1)=1.0;%0.9;
        else factormax(ind1)=1.0;%1.1;
        end
    end
    if nargin==5
        hndl1=varargin{5};
    else
        hndl1=gca;
    end
    Nticks_x=5;
    Nticks_y=5;
    color_left='red';
    color_right='blue';
    linestyle{1}='-';
    linestyle{2}='--';
    [AX, H1, H2] = plotyy(hndl1,x,y1,x,y2);
%     set(AX(1),'xlim',[min(x) max(x)],'Xtick',roundP([min(x):(max(x)-min(x))/(Nticks_x-1):max(x)],2));
%     set(AX(2),'xlim',[min(x) max(x)],'Xtick',roundP([min(x):(max(x)-min(x))/(Nticks_x-1):max(x)],2));
    set(AX(1),'xlim',[min(x) max(x)],'Xtick',DecadeTicks([min(x) max(x)]));
    set(AX(2),'xlim',[min(x) max(x)],'Xtick',DecadeTicks([min(x) max(x)]));
    if min(min(y1))==0 && max(max(y1))==0
        set (AX(1), 'Ycolor', color_left,'xgrid','on','Ylim',[-1 1],'Ytick',roundP([-1:2/(Nticks_y-1):1],2))
    else
        set (AX(1), 'Ycolor', color_left,'xgrid','on','Ylim',[min(min(y1))*factormin(1),max(max(y1))*factormax(1)],'Ytick',roundP([min(min(y1))*factormin(1):(max(max(y1))*factormax(1)-min(min(y1,1))*factormin(1))/(Nticks_y-1):max(max(y1))*factormax(1)],2));
    end
    if min(min(y2))==0 && max(max(y2))==0
        set (AX(2), 'Ycolor', color_right,'Ylim',[-1 1],'Ytick',roundP(-1:2/(Nticks_y-1):1,2))
    else
        set (AX(2), 'Ycolor', color_right,'Ylim',[min(min(y2))*factormin(2),max(max(y2))*factormax(2)],'Ytick',roundP([min(min(y2))*factormin(2):(max(max(y2))*factormax(2)-min(min(y2))*factormin(2))/(Nticks_y-1):max(max(y2))*factormax(2)],2));
    end
    for ind1=1:length(H1)
        set(H1(ind1), 'Linestyle', linestyle{ind1}, 'Linewidth', 2, 'Color', color_left);
    end
    for ind1=1:length(H2)
        set(H2(ind1), 'Linestyle', linestyle{ind1}, 'Linewidth', 2, 'Color', color_right);
    end
    % set(AX,'xlim',[nu_lim1, nu_lim2]);
    %xlim ([nu_lim1 nu_lim2]);
    setaxesP(AX)
    if nargin>=4
        texts=varargin{4};
        xlabel(AX(1),texts{1}, 'fontname', 'verdana', 'fontsize', 14);
        set (get(AX(1), 'Ylabel'), 'String', texts{2}, 'fontname', 'verdana', 'fontsize', 14, 'color', color_left);
        set (get(AX(2), 'Ylabel'), 'String', texts{3}, 'fontname', 'verdana', 'fontsize', 14, 'color', color_right);
        set (get(AX(1), 'Title'), 'String', texts{4}, 'fontname', 'verdana', 'fontsize', 10, 'color', 'k');
    end
end
end
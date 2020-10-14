function varargout = plotyyP(input)

% input is a struct, with the following mandatory fields:
% x,y1,y2,y1lim,y2lim
% optional fields:
% xlim,xlabel,y1_label,y2_label,xtickpos,Nticks_x,Nticks_y,color_left,color_right,lwidth,fntsz,fnt,Nround

% 'x' has to be of the size Nx1, and y1 has to be NxM, y2 has to be NxK.

x=input.x;
Ny1 = size(input.y1,2);
Ny2 = size(input.y2,2);
y(:,1:Ny1)=input.y1;
y(:,Ny1+1:Ny1+Ny2)=input.y2;
for ind1=1:Ny1+Ny2
    if min(y(:,ind1))<0;
        factormin(ind1)=1.0;%1.1;
    else factormin(ind1)=1.0;%0.9;
    end
    if max(y(:,ind1))<0;
        factormax(ind1)=1.0;%0.9;
    else factormax(ind1)=1.0;%1.1;
    end
end
if isfield(input,'hax')
    hndl1=input.hax;
else
    hndl1=gca;
end
if isfield(input,'Nticks_x')
    Nticks_x=input.Nticks_x;
else
    Nticks_x=5;
end
if isfield(input,'Nticks_y')
    Nticks_y=input.Nticks_y;
else
    Nticks_y=5;
end
if isfield(input,'color_left')
    color_left=input.color_left;
else
    color_left='black';
end
if isfield(input,'color_right')
    color_right=input.color_right;
else
    color_right='blue';
end
if isfield(input,'lwidth')
    lwidth=input.lwidth;
else
    lwidth=2;
end
if isfield(input,'fnt')
    fnt=input.fnt;
else
    fnt='Verdana';
end
if isfield(input,'fntsz')
    fntsz=input.fntsz;
else
    fntsz=14;
end
if isfield(input,'Nround')
    Nround=input.Nround;
else
    Nround=2;
end
[AX, H1, H2] = plotyy(hndl1,x,y(:,1:Ny1),x,y(:,Ny1+1:Ny1+Ny2));
varargout{1} = [AX, H1, H2];
set (AX(1), 'Ycolor', color_left,'xgrid','off','Ylim',input.y1lim,'Ytick',roundP([input.y1lim(1):(input.y1lim(2)-input.y1lim(1))/(Nticks_y-1):input.y1lim(2)],Nround));
set (AX(2), 'Ycolor', color_left,'xgrid','off','Ylim',input.y2lim,'Ytick',roundP([input.y2lim(1):(input.y2lim(2)-input.y2lim(1))/(Nticks_y-1):input.y2lim(2)],Nround));

% if min(y(:,1))==0 && max(y(:,2)==0)
%     set (AX(1), 'Ycolor', color_left,'xgrid','off','Ylim',[-1 1],'Ytick',roundP([-1:2/(Nticks_y-1):1],2))
% else
%     set (AX(1), 'Ycolor', color_left,'xgrid','off','Ylim',[min(y(:,1))*factormin(1),max(y(:,1))*factormax(1)],'Ytick',roundP([min(y(:,1))*factormin(1):(max(y(:,1))*factormax(1)-min(y(:,1))*factormin(1))/(Nticks_y-1):max(y(:,1))*factormax(1)],2));
% end
% if min(y(:,2))==0 && max(y(:,2)==0)
%     set (AX(2), 'Ycolor', color_right,'Ylim',[-1 1],'Ytick',roundP([-1:2/(Nticks_y-1):1],2))
% else
%     set (AX(2), 'Ycolor', color_right,'Ylim',[min(y(:,2))*factormin(2),max(y(:,2))*factormax(2)],'Ytick',roundP([min(y(:,2))*factormin(2):(max(y(:,2))*factormax(2)-min(y(:,2))*factormin(2))/(Nticks_y-1):max(y(:,2))*factormax(2)],2));
% end
set(H1, 'Linestyle', '-', 'Linewidth', lwidth, 'Color', color_left);
set(H2, 'Linestyle', '-', 'Linewidth', lwidth, 'Color', color_right);
set(AX(1),'ycolor',color_left);
set(AX(2),'ycolor',color_right);
setfigP(AX(1),[fntsz lwidth]);
setfigP(AX(2),[fntsz lwidth])
if isfield(input,'xlabel')
    xlabel(AX(1),input.xlabel,'fontname',fnt);
end
if isfield(input,'xtickpos')
    set(AX(1),'Xtick',input.xtickpos);
    set(AX(2),'Xtick',input.xtickpos);
end
if isfield(input,'xlim')
    set(AX(1),'xlim',input.xlim,'Xtick',roundP([input.xlim(1):(input.xlim(2)-input.xlim(1))/(Nticks_x-1):input.xlim(2)],2));
    set(AX(2),'xlim',input.xlim,'Xtick',roundP([input.xlim(1):(input.xlim(2)-input.xlim(1))/(Nticks_x-1):input.xlim(2)],2));
end
if isfield(input,'y1_label')
    set(get(AX(1),'Ylabel'), 'String',input.y1_label,'fontname',fnt,'color',color_left);
end
if isfield(input,'y2_label')
    set(get(AX(2),'Ylabel'), 'String',input.y2_label,'fontname',fnt,'color',color_right);
end
if isfield(input,'title')
    set(get(AX(1),'Title'), 'String',input.title,'fontname',fnt,'color','k');
end
end
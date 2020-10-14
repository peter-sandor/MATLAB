function plotsave2D(frq3,frq1,data2D,str_title,varargin)

cont_max=max(max(real(data2D)));
% cont_min=mean(mean(real(data_2D(1:50,1:50))))+cont_max/25;
cont_min=min(min(real(data2D)));
cont_N=20;
hndl3=figure;contour(frq3,frq1,real(data2D),cont_min:(cont_max-cont_min)/(cont_N-1):cont_max)
xlabel('\nu_3 [PHz]')
ylabel('\nu_1 [PHz]')
% xlim([min(frq3) max(frq3)])
% ylim([min(frq1) max(frq1)])
grid on
if nargin>=5
    if isnumeric(varargin{1})
        xlim(varargin{1})
        ylim(varargin{1})
    elseif nargin==6
        xlim(varargin{2}(1,:))
        ylim(varargin{2}(2,:))
    else    xlim([1.135 1.170])
            ylim([1.135 1.170])
    end
    
end
title(str_title);
colorbar
if nargin>=5
    if isstr(varargin{1})
        saveas(hndl3,varargin{1});
    end
end
close(hndl3)
end
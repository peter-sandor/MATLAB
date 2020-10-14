function coords = get_data_coords(varargin)
if nargin==1
    hax=gca;
    value=varargin{1};
elseif nargin==2
    hax=varagin{1};
    value=varargin{2};
end
% xax_data=get(get(hax,'Children'),'xdata');
ax_pos=get(hax,'position');
xax_xlim=get(hax,'xlim');
xax_scale=get(hax,'xscale');
yax_xlim=get(hax,'ylim');
yax_scale=get(hax,'yscale');
if strcmp(xax_scale,'linear')
    xcoord=ax_pos(3)*(value(1)-xax_xlim(1))/(xax_xlim(2)-xax_xlim(1))+ax_pos(1);
elseif strcmp(xax_scale,'log')
    xcoord=ax_pos(3)*(log10(value(1))-log10(xax_xlim(1)))/(log10(xax_xlim(2))-log10(xax_xlim(1)))+ax_pos(1);
end
if strcmp(yax_scale,'linear')
    ycoord=ax_pos(4)*(value(2)-yax_xlim(1))/(yax_xlim(2)-yax_xlim(1))+ax_pos(2);
elseif strcmp(yax_scale,'log')
    ycoord=ax_pos(4)*(log10(value(2))-log10(yax_xlim(1)))/(log10(yax_xlim(2))-log10(yax_xlim(1)))+ax_pos(2);
end
coords=[xcoord ycoord];
end
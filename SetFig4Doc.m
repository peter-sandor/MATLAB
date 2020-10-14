function SetFig4Doc(hfig,input)

width=input(1); % figure width in inches
fntsz=input(2); % fontsize in points
lnthck=input(3); % line thickness in points
SetFigWidth(width);
hAllAxes = findobj(hfig,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg); % All axes which are not legends

for ind1=1:length(hAxes)
    setfigP(hAxes(ind1),[fntsz lnthck]);
end
end
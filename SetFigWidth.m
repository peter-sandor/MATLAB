function SetFigWidth(width)

hfig1=gcf;
hfig1.Units = 'inches';
aspect_ratio=hfig1.Position(4)/hfig1.Position(3); % height/width
hfig1.Position(3)=width;
hfig1.Position(4)=width*aspect_ratio;
end
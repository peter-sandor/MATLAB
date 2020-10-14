function print2eps(varargin)
if nargin==0
    hfig=gcf;
else
    hfig=varargin{1};
end

hfig.PaperUnits = 'inches';
hfig.PaperPositionMode = 'auto';
hfig.PaperPosition=hfig.Position;
filename=uigetfile('*.eps');
if filename~=0
    print(hfig, '-opengl', '-depsc', '-r300',filename);
    disp('Print successful.')
end
end
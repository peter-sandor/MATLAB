function [ax,hl1,hl2] = plotxxP(x1,x2,y,xlabels,ylabels);
%PLOTXX - Create graphs with x axes on both top and bottom 
%
%Similar to PLOTYY, but ...
%the independent variable is on the y-axis, 
%and both dependent variables are on the x-axis.
%
%Syntax: [ax,hl1,hl2] = plotxx(x1,y,x2,y2,xlabels,ylabels);
%
%Inputs:  X1,Y1 are the data for the first line (black)
%         X2,Y2 are the data for the second line (red)
%         XLABELS is a cell array containing the two x-labels
%         YLABELS is a cell array containing the two y-labels
%
%The optional output handle graphics objects AX,HL1,HL2
%allow the user to easily change the properties of the plot.
%
%Example: Plot temperature T and salinity S 
%         as a function of depth D in the ocean
%
%D = linspace(-100,0,50);
%S = linspace(34,32,50);
%T = 10*exp(D/40);
%xlabels{1} = 'Temperature (C)';
%xlabels{2} = 'Salinity';
%ylabels{1} = 'Depth(m)';
%ylabels{2} = 'Depth(m)';
%[ax,hlT,hlS] = plotxx(T,D,S,D,xlabels,ylabels);


%The code is inspired from page 10-26 (Multiaxis axes)
%of the manual USING MATLAB GRAPHICS, version 5.
%
%Tested with Matlab 5.3.1 and above on PCWIN

%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%November 1997; Last revision: 01-Nov-2001

if x1(end)-x1(1)>0
    ax1_dir = 'normal';
else
    ax1_dir = 'reverse';
end

if x2(end)-x2(1)>0
    ax2_dir = 'normal';
else
    ax2_dir = 'reverse';
end

if nargin < 4
   error('Not enough input arguments')
elseif nargin==4
   %Use empty strings for the xlabels
   xlabels{1}=' '; xlabels{2}=' '; ylabels{1}=' '; ylabels{2}=' ';
elseif nargin==5
   %Use empty strings for the ylabel
   ylabels{1}=' '; ylabels{2}=' ';
elseif nargin > 6
   error('Too many input arguments')
end

if length(ylabels) == 1
   ylabels{2} = ' ';
end

if ~iscellstr(xlabels) 
   error('Input xlabels must be a cell array')
elseif ~iscellstr(ylabels) 
   error('Input ylabels must be a cell array')
end

val_pad = 20;
ax_pos = [0.12 0.15 0.85 0.65];
ylim_vals = [min(y) + (min(y)-max(y))/val_pad, max(y) + (max(y)-min(y))/val_pad];

ax(1) = axes('Position',ax_pos,...
    'xdir',ax1_dir,...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k');

hl1 = line(x1,y,'Parent',ax(1),'linestyle','none','Marker','none','Markersize',10,'Color','k');

ax(2) = axes('Position',ax_pos,...
    'xdir',ax1_dir,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','w',...
    'XColor','r','YColor','k',...
    'XTick',x1,...
    'XTickLabels',roundP(x2,2));

hl2 = line(x1,y,'Parent',ax(2),'linestyle','none','Marker','o','Markersize',10,'Color','k');

xlim(ax(1),[min(x1) + (min(x1)-max(x1))/val_pad, max(x1) + (max(x1)-min(x1))/val_pad]);
xlim(ax(2),[min(x1) + (min(x1)-max(x1))/val_pad, max(x1) + (max(x1)-min(x1))/val_pad]);
ylim(ax(1),ylim_vals);
% setlines;
set(ax(2),'YTick',[]);

set(ax,'box','off')
linkaxes(ax,'y');
setfigP(ax(1));
setfigP(ax(2));
%label the two x-axes
set(get(ax(1),'xlabel'),'string',xlabels{1})
set(get(ax(2),'xlabel'),'string',xlabels{2})
set(get(ax(1),'ylabel'),'string',ylabels{1})
set(get(ax(2),'ylabel'),'string',ylabels{2})

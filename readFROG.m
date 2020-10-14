function readFROG

[FileName,PathName] = uigetfile({'*.?fa','FROG data & axes';'*.*','All Files (*.*)';'*.ifa','FROG data';'*.xfa','FROG axes'},'Pick a FROG data and a FROG axes file','MultiSelect', 'on');
flagi=-1;
flagx=-1;
for ind1=1:length(FileName)
    if isempty(strfind(FileName{ind1},'.ifa'))==0
        flagi=ind1;
    elseif isempty(strfind(FileName{ind1},'.xfa'))==0
        flagx=ind1;
    end
end

if (flagi>0)
%     data=load(strcat(PathName,FileName{flagi}));
    [truncated,header]=deletelines(FileName{flagi},5,'FROG.txt');
    S=importdata(strcat(PathName,truncated));
    FROGdata=dlmread('FROG.txt');;
end

if (flagx>0)
    temp=load(strcat(PathName,FileName{flagx}));
end

if length(unique(temp(:,2)))==length(temp(:,2))
    axes{1}=temp(1:header(1),1);
    axes{2}=temp(:,2);
end

hndl1=figure;imagesc(FROGdata)
title('Select region of interest')
temp=ginput(2);
index1_start=uint16(round(min(temp(:,2))));
index1_end=uint16(round(max(temp(:,2))));
index2_start=uint16(round(min(temp(:,1))));
index2_end=uint16(round(max(temp(:,1))));
if index1_start < 1
    index1_start=1;
end
if index2_start < 1
    index2_start=1;
end
if index1_end > size(FROGdata,1)
    index1_end=size(FROGdata,1);
end
if index2_end > size(FROGdata,2)
    index2_end=size(FROGdata,2);
end 
FROGdata=FROGdata(index1_start:index1_end,index2_start:index2_end).';
axes{1}=axes{1}(index1_start:index1_end);
axes{2}=axes{2}(index2_start:index2_end);
% close(hndl1)
contour(axes{1},axes{2},FROGdata,22)
ylabel('Wavelength [nm]')
xlabel('Delay [fs]')
title(strcat(PathName,FileName{flagi}),'Interpreter','None')
colorbar

saveas(hndl1,strcat(PathName,FileName{flagi}(1:length(FileName{flagi})-4),'.fig'))
end
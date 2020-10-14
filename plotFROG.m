function plotFROG(varargin)

if nargin==1
    origindir=[pwd '\'];
    PathName=varargin{1};
    cd(PathName);
    files=dir('*.ifa');
    if ~isempty(files)
        FileName=files(1).name;
    end
else
%     [FileName,PathName] = uigetfile({'*.*','All Files (*.*)';'*.ifa','FROG data'},'Pick a FROG data file (.ifa)','MultiSelect', 'off');
    files=dir('*.ifa');
    PathName='.\';
    if ~isempty(files)
        FileName=files(1).name;
    end
end

if ~isempty(FileName)
    [truncated,header]=deletelines(FileName,5,'FROG.txt');
    S=importdata(strcat(PathName,truncated));
    FROGdata=dlmread('FROG.txt').';

    Ndelay=header(1);
    delay_step=header(3); % [fs/step]
    delay_min=-(Ndelay/2-1)*delay_step;
    delay_max=Ndelay/2*delay_step;
    delay=delay_min:delay_step:delay_max;

    Nlambda=header(2);
    lambda_step=header(4); % [nm/step]
    lambda_center=header(5); % [nm]
    lambda_min=lambda_center-lambda_step*Nlambda/2;
    lambda_max=lambda_min+(Nlambda-1)*lambda_step;
    lambda=lambda_min:lambda_step:lambda_max;

    axes{1}=delay;
    axes{2}=lambda;

    % hndl1=figure;
    % imagesc(FROGdata)
    % title('Select region for background subtraction')
    % temp=ginput(2);
    % index1_start=uint16(round(min(temp(:,2))));
    % index1_end=uint16(round(max(temp(:,2))));
    % index2_start=uint16(round(min(temp(:,1))));
    % index2_end=uint16(round(max(temp(:,1))));
    % if index1_start < 1
    %     index1_start=1;
    % end
    % if index2_start < 1
    %     index2_start=1;
    % end
    % if index1_end > size(FROGdata,1)
    %     index1_end=size(FROGdata,1);
    % end
    % if index2_end > size(FROGdata,2)
    %     index2_end=size(FROGdata,2);
    % end 
    index1_start=1;
    index1_end=20;
    index2_start=1;
    index2_end=20;
    meanvalue=mean(mean(FROGdata(index1_start:index1_end,index2_start:index2_end)));

    % title('Select region of interest')
    % temp=ginput(2);
    % index1_start=uint16(round(min(temp(:,2))));
    % index1_end=uint16(round(max(temp(:,2))));
    % index2_start=uint16(round(min(temp(:,1))));
    % index2_end=uint16(round(max(temp(:,1))));
    % if index1_start < 1
    %     index1_start=1;
    % end
    % if index2_start < 1
    %     index2_start=1;
    % end
    % if index1_end > size(FROGdata,1)
    %     index1_end=size(FROGdata,1);
    % end
    % if index2_end > size(FROGdata,2)
    %     index2_end=size(FROGdata,2);
    % end
    [X,Y]=meshgrid(1:size(FROGdata,2),1:size(FROGdata,1));
    temp=(FROGdata>2*meanvalue);
    index1_start=min(Y(temp));
    index1_end=max(Y(temp));
    index2_start=min(X(temp));
    index2_end=max(X(temp));

    FROGdata=(FROGdata(index1_start:index1_end,index2_start:index2_end)-meanvalue).'/max(max(FROGdata(index1_start:index1_end,index2_start:index2_end)-meanvalue));
    axes{1}=axes{1}(index1_start:index1_end);
    axes{2}=axes{2}(index2_start:index2_end);
    % close(hndl1)

    hndl2=figure;
    hax(1)=subplot(223);
    contour(axes{1},axes{2},FROGdata,22)
    set(hax(1),'yaxislocation','right')
    xlim([min(axes{1}) max(axes{1})]);
    ylim([min(axes{2}) max(axes{2})]);
    ylabel('Wavelength [nm]')
    xlabel('Delay [fs]')
    title(strcat(PathName,FileName),'Interpreter','None')
    colorbar('West')

    subplot(221);
    plot(axes{1},sum(FROGdata,1),'k','LineWidth',2)
    xlim([min(axes{1}) max(axes{1})]);
    ylim([0 max(sum(FROGdata,1))])
    title('Autocorrelation')

    subplot(224)
    plot(sum(FROGdata,2),axes{2},'k','LineWidth',2)
    ylim([min(axes{2}) max(axes{2})]);
    xlim([0 max(sum(FROGdata,2))])
    title('Spectrum')

    saveas(hndl2,strcat(PathName,FileName(1:length(FileName)-4),'.fig'));
    close(hndl2)

end
if nargin==1
    cd(origindir);
end
end
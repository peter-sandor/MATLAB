function [path_file, totalint] = procbin(h)
        [FileName,PathName] = uigetfile;
        path_file=strcat(PathName,FileName);
        if isempty(findstr(FileName,'.bin'))==0
            data=loadbinaryimage(path_file);
        elseif isempty(findstr(FileName,'.jpg'))==0
            data=double(imread(path_file));
        else data=[];
            totalint=NaN;
            return
        end
        % construct subset of image that contains measured data
        center=[527;399]; % x,y
        radius=350;
        index_fig=zeros(size(data));
        for ind1=1:size(data,2)
            for ind2=1:size(data,1)
                if ((ind1-center(1))^2+(ind2-center(2))^2)<=radius^2
                    index_fig(ind2,ind1)=1;
                end
            end
        end
        index_fig=logical(index_fig);
        % integrate over the subset of image, and place image in the
        % appropriate axes
        totalint=sum(sum(data(index_fig)));
        axes(h);
        data_plot=zeros(size(data));
        data_plot(index_fig)=data(index_fig);
        imagesc(data_plot);
        colorbar
end
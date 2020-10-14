function varargout = centroid(varargin)
% This code takes an image (e.g. acquired with VMI machine) and determines
% the locations of electron (or ion) hits, and returns with the coordinates.
% The idea is to first threshold the image, and then sort the pixels with
% nonzero values into lists (or sets) based on topological connectedness.
% Each list will consist of pixels that are supposed to be topologically
% connected, and hence belong to the same hit. Once the lists for the
% different hits are established, the coordinate of a specific hit can be calculated by
% taking the intensity-averaged mean of the pixel coordinates on its
% corresponding list.
% Code execution time scales as the square of the number of pixels to be
% sorted, so for many hits  on the image, expect long runtime.

% usage: coords_out = centroid
%                   in this case an image is read in from a file and a
%                   default threshold value is used. Code returns with an
%                   Nx2 array, where in each row coordinates of a hit are
%                   stored in the format [y x] = [column row]
%        ... = centroid(img,thrs)
%                   here, 'img' is either a 2D array containing the image,
%                   or a filename (with full path) to the file.
%                   'thrs' is the threshold value to be used

radius=2; % distance based on which it is determined whether two pixels are topologically connected or not (are they close enough to each other)
minlistlength=4; % this number sets the threshold for a check later on in the code that tells whether each list contains a number of pixels that is reasonable large
filtersize=60;
if nargin==0
    [FileName,PathName] = uigetfile({'*.jpg','JPEG image';'*.png','PNG image';'*.bmp','BMP image';'*.*','All Files (*.*)'},'Pick an image file');
    path_file=strcat(PathName,FileName);
    pic=monopic(imread(path_file));
    thrs=10;
elseif nargin==1
    if ischar(varargin{1})
        path_file=varargin{1};
        pic=monopic(imread(path_file));
    elseif isnumeric(varargin{1})
        pic=varargin{1};
    else
        disp('Invalid input.')
        varargout{1}=[];
        varargout{2}=[];
        return
    end
    thrs=10;
elseif nargin==2
    if ischar(varargin{1})
        path_file=varargin{1};
        pic=monopic(imread(path_file));
    elseif isnumeric(varargin{1})
        pic=varargin{1};
    else
        disp('Invalid input.')
        varargout{1}=[];
        varargout{2}=[];
        return
    end
    thrs=varargin{2};
end
% N=size(pic);

pts=thrs_img(pic,thrs,radius,minlistlength,filtersize); % threshold image and enumerate the points surviving the process --> 'pts'

if size(pts,1)>1    % if the list is not empty, we start sorting its elements into lists, which will corresponds to the different hits on the detector
%     pic=double(pic);
    M=size(pts,1);
    % Now create connection matrix with size M x M, logical type. Two pixels are connected
    % if their distance is smaller than the value of 'radius'.
    ind2=1;
    ind_keep=[];
    for ind1=1:M
        temp=((pts(ind1,1)-pts(:,1)).^2+(pts(ind1,2)-pts(:,2)).^2)<=radius^2;
        if sum(temp)>minlistlength % filter elements of the list additionally for standalone pixels with no neighbours that survived the thresholding. These are surely not valid hits.
            Mconnect(:,ind2)=logical(temp);
            ind2=ind2+1;
            ind_keep=[ind_keep; ind1];
        end
    end
        
%     %filter out "junk" single pixels
%     index1=(sum(double(Mconnect),1)<=minlistlength);
%     pts(index1,:)=[];
%     Mconnect(index1,:)=[];
%     Mconnect(:,map2rowvec(index1))=[];

    if ~isempty(Mconnect)
        Mconnect=Mconnect(ind_keep,:);
        pts=pts(ind_keep,:);
        % The next part uses the connection matrix to sort the pixels.
        for ind1=1:size(Mconnect,1)-1
            if ind1<=size(Mconnect,1)/2
                range=min(ind1-1,50);
            else
                range=min(size(Mconnect,1)-ind1,50);
            end
            summed(ind1)=sum(sum(double(Mconnect(ind1-range:ind1,ind1+1:ind1+range))));
        %     disp([num2str(ind1) ' ' num2str(range)]);
        end

        sumlogic=true([length(summed) 1]);
        sumlogic(summed<=5)=false;

        flag=false;
        ind2=1;
        for ind1=1:length(sumlogic)
            if sumlogic(ind1) && ~flag % start of range of pixels for a hit
                flag=true;
                indices(ind2,1)=ind1;
            elseif ~sumlogic(ind1) && flag % end of range of pixels
                indices(ind2,2)=ind1;
                flag=false;
                ind2=ind2+1;
            end
        end

        ind2=1;
        for ind1=1:size(indices,1)
            S0{ind1}=pts(indices(ind1,1):indices(ind1,2),:);
            M0{ind1}=Mconnect(indices(ind1,1):indices(ind1,2),indices(ind1,1):indices(ind1,2));

            while ~isempty(M0{ind1})
            % find first-order connections
            member0=logical(M0{ind1}(:,1)); % index vector (column)
            % find second-order connections
            member1=any(M0{ind1}(member0,:),1);
            member1=map2colvec(member1); % index vector (column)
            member_diff=xor(member0,member1);
            M0{ind1}(member0,:)=false;
            M0{ind1}(:,map2rowvec(member0))=false;
            member_tot=member1;
            % check and follow up on higher-order connections
            while any(member_diff)
                member0=member1;
                member1=any(M0{ind1}(member_diff,:),1);
                member1=map2colvec(logical(member1)); % index vector (column)
                member_diff=xor(member0,member1);
                member_tot=or(member_tot,member1);
                M0{ind1}(member0,:)=false;
                M0{ind1}(:,map2rowvec(member0))=false;
            end
            S{ind2}=S0{ind1}(member_tot,:); % 'member_tot' contains the indices for all the connected pixels, here create a separate list for them. 
            S0{ind1}(member_tot,:)=[]; % Delete elements that have already been taken into account. This speeds up the sorting.
            M0{ind1}(member_tot,:)=[];
            M0{ind1}(:,map2rowvec(member_tot))=[];
            ind2=ind2+1;
            end

        end
        % Done with creating the lists, next step is to see if they are 'real'
        % hits. In case yes, their center coordinates are calculated.
        hits=0;
        cntr=[];
        for ind1=1:length(S)
            if length(S{ind1})>minlistlength % checks whether each list contains a number of pixels that is reasonable large
                                             % this helps filter out single pixel "defects" with large illumination value on them, which may result in false hits
                hits=hits+1;
                cntr(hits,:)=[0 0];
                for ind2=1:length(S{ind1})
                    cntr(hits,:)=cntr(hits,:)+S{ind1}(ind2,3)*S{ind1}(ind2,1:2);
                end
                cntr(hits,:)=cntr(hits,:)/sum(S{ind1}(:,3));
            end
        end
        
    else
        cntr=[];
    end
else
    cntr=[];
end
varargout{1}=cntr; % return with a list of coordinates
end
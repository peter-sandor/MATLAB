function varargout = centroid_old(varargin)
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
%                   stored in the format [x y] = [row column]
%        ... = centroid(img,thrs)
%                   here, 'img' is either a 2D array containing the image,
%                   or a filename (with full path) to the file.
%                   'thrs' is the threshold value to be used


radius=2; % distance based on which it is determined whether two pixels are topologically connected or not (are they close enough to each other)
minlistlength=20; % this number sets the threshold for a check later on in the code that tells whether each list contains a number of pixels that is reasonable large
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
N=size(pic);
pic(pic<=thrs)=0; % thresholding is done here
if max(max(pic))>0 % if there are pixels with nonzero values still, we start sorting them into the lists
    pic=double(pic);
    [y,x]=find(pic);
    pts(:,1)=x;
    pts(:,2)=y;
    M=size(pts,1);
    S{1}=pts(1,:); % create the first list with one element
    for ind1=1:M
        member1=0; 
        for ind3=1:length(S)
            if logical(max(findvec(S{ind3},pts(ind1,:)))) % check whether point #1 is on any list
                member1=ind3;
            end
        end
        if member1==0 % if not, create new list, and make point #1 to be its first and only element so far
            member1=length(S)+1;
            S{member1}=pts(ind1,:);
        end
        if length(S)==3
            ind1;
            ind2;
        end
        for ind2=(ind1+1):M
            if (sqrt((pts(ind1,1)-pts(ind2,1))^2+(pts(ind1,2)-pts(ind2,2))^2)<=radius) % check whether point #1 and #2 belong to the same set
                member2=0;
                for ind3=1:length(S)
                    if logical(max(findvec(S{ind3},pts(ind2,:)))); % point #1 and point #2 belong together; now we look if point #2 is an element of any list already
                        member2=ind3;
                    end
                end
                if member2==0 % if point #2 is not listed, put it on the same as #1
                    S{member1}=vertcat(S{member1},pts(ind2,:));
                elseif member2~=0 && (member1~=member2) % if it is listed, merge the lists of #1 and #2
                    S{member1}=vertcat(S{member1},S{member2});
                    S(member2)=[];
                    if member1>member2
                        member1=member1-1;
                    end
                end
            end
        end
    end
    hits=0;
    cntr=[];
    for ind1=1:length(S)
        if length(S{ind1})>minlistlength % checks whether each list contains a number of pixels that is reasonable large
                                         % this helps filter out single pixel "defects" with large illumination value on them, which may result in false hits
            hits=hits+1;
            cntr(hits,:)=[0 0];
            weights=zeros([length(S{ind1}) 1]);
            for ind2=1:length(S{ind1})
                weights(ind2)=pic(S{ind1}(ind2,2),S{ind1}(ind2,1));
                cntr(hits,:)=cntr(hits,:)+weights(ind2)*S{ind1}(ind2,:);
            end
            cntr(hits,:)=cntr(hits,:)/sum(weights);
        end
    end
else cntr=[];
    hits=0;
end
% varargout{2}=hits;
varargout{1}=cntr;
% disp(num2str(toc));
end
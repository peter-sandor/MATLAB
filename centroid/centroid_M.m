function varargout = centroid_M(varargin)
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
%        ... = centroid(img,params)
%                   here, 'img' is either a 2D array containing the image,
%                   or a filename (with full path) to the file.
%                   params is a struct with the following fields:
%                   'thrs': the threshold value to be used
%                   'radius': distance based on which it is determined whether two pixels are topologically connected or not (are they close enough to each other)
%                   'minlistlength': this number sets the threshold for a check later on in the code that tells whether each list contains a number of pixels that is reasonable large
%                   'filtersize': how many neighbors should the initial filtering routine check?
%                   'max_coord_diff': this number controls a filtering condition, a good value is 25

if nargin==0
    [FileName,PathName] = uigetfile({'*.jpg','JPEG image';'*.png','PNG image';'*.bmp','BMP image';'*.*','All Files (*.*)'},'Pick an image file');
    path_file=strcat(PathName,FileName);
    pic=monopic(imread(path_file));
    thrs=10;
    radius=2;
    minlistlength=4;
    filtersize=60;
    swap_xy = 0;
else
    pic=varargin{1};
    thrs=varargin{2}.thrs;
    radius=varargin{2}.radius; 
    minlistlength=varargin{2}.minlistlength;
    filtersize=varargin{2}.filtersize;
    max_coord_diff=varargin{2}.max_coord_diff;
    swap_xy=varargin{2}.swap_xy;
end

pts=thrs_img(pic,thrs,radius,minlistlength,filtersize,swap_xy); % threshold image and enumerate the points surviving the process --> 'pts'
if size(pts,1)>1    % if the list is not empty, we start sorting its elements into lists, which will corresponds to the different hits on the detector
    Mconnect=create_Mconnect(pts,radius,minlistlength);
    indices=presort_pixels(Mconnect,filtersize);
    
        S=[];
        for ind1=1:size(indices,1)
            S0=pts(indices(ind1,1):indices(ind1,2),:);
            x_diff=max(S0(:,1))-min(S0(:,1));
            y_diff=max(S0(:,2))-min(S0(:,2));
            if ~isempty(indices)
                M0=Mconnect(indices(ind1,1):indices(ind1,2),indices(ind1,1):indices(ind1,2));
                if (max(x_diff,y_diff)>max_coord_diff)
                    temp = cluster_pixels_M(S0,M0);
                    S = cat(2,S,temp);
                else
                    temp2{1}=S0;
                    S = cat(2,S,temp2);
                end
            end
        end
        % Done with creating the lists, next step is to see if they are 'real'
        % hits. In case yes, their center coordinates are calculated.
        coords = calculate_coordinates(S,minlistlength);
else
    coords=[];
end
varargout{1}=coords; % return with a list of coordinates
varargout{2}=S;
end

function varargout = thrs_img(img,thrs,radius,minlistlength,filtersize,swap_xy)
% Create a list of coordinates which are above the given threshold ('thrs')

N=size(img);
ind3=1;
pts=zeros([10000 3]);
for ind1=2:N(1)-1
    for ind2=2:N(2)-1
%         if (img(ind1,ind2)+img(ind1-1,ind2)+img(ind1,ind2-1)+img(ind1+1,ind2)+img(ind1,ind2+1))>5*thrs
         if img(ind1,ind2)>thrs
            pts(ind3,:)=[ind1 ind2 img(ind1,ind2)];
            ind3=ind3+1;
        end
    end
end
pts(pts(:,1)==0,:)=[];
if size(pts,1)>1
    pts2=[];
    ind2=1;
    Npts=size(pts,1);
    for ind1=1:Npts
        range_min=max(1,ind1-filtersize);
        range_max=min(Npts,ind1+filtersize);
        temp=double(((pts(ind1,1)-pts(range_min:range_max,1)).^2+(pts(ind1,2)-pts(range_min:range_max,2)).^2)<=radius^2);
        if sum(temp)>minlistlength
            pts2(ind2,:)=pts(ind1,:);
            ind2=ind2+1;
        end
    end
else
    pts2=[];
end
if swap_xy
    varargout{1}=sortrows([pts2(:,2) pts2(:,1) pts2(:,3)],2);
else    
    varargout{1}=pts2;
end
end

function Mconnect = create_Mconnect(pts,radius,filtersize)
% Create connection matrix with size M x M, logical type. Two pixels are connected
% if their distance is smaller than the value of 'radius'.
M=size(pts,1);
% ind2=1;
% ind_keep=[];
Mconnect=zeros([M M]);
for ind1=1:M
    temp=((pts(ind1,1)-pts(:,1)).^2+(pts(ind1,2)-pts(:,2)).^2)<=radius^2;
%     if sum(temp)>minlistlength % filter elements of the list additionally for standalone pixels with no neighbours that survived the thresholding. These are surely not valid hits.
        Mconnect(:,ind1)=logical(temp);
%         ind2=ind2+1;
%         ind_keep=[ind_keep; ind1];
%     end
end
%         Mconnect=Mconnect(ind_keep,:);
%         pts=pts(ind_keep,:);
end

function indices = presort_pixels(Mconnect,filtersize)

summed=zeros([size(Mconnect,1)-1 1]);
for ind1=1:size(Mconnect,1)-1
    if ind1<=size(Mconnect,1)/2
        range=min(ind1-1,filtersize);
    else
        range=min(size(Mconnect,1)-ind1,filtersize);
    end
    summed(ind1)=sum(sum(double(Mconnect(ind1-range:ind1,ind1+1:ind1+range))));
%     disp([num2str(ind1) ' ' num2str(range)]);
end

sumlogic=true([length(summed) 1]);
sumlogic(summed<=5)=false;

% indices=zeros([length(sumlogic) 2]);
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
if ~exist('indices','var');
    indices=[];
end
end

function S = cluster_pixels_M(pts,Mconnect)

ind1=1;
while ~isempty(Mconnect)
    % find first-order connections
    member0=logical(Mconnect(:,1)); % index vector (column)
    % find second-order connections
    member1=any(Mconnect(member0,:),1);
    member1=map2colvec(member1); % index vector (column)
    member_diff=xor(member0,member1);
    Mconnect(member0,:)=false;
    Mconnect(:,map2rowvec(member0))=false;
    member_tot=member1;
    % check and follow up on higher-order connections
    while any(member_diff)
        member0=member1;
        member1=any(Mconnect(member_diff,:),1);
        member1=map2colvec(logical(member1)); % index vector (column)
        member_diff=xor(member0,member1);
        member_tot=or(member_tot,member1);
        Mconnect(member0,:)=false;
        Mconnect(:,map2rowvec(member0))=false;
    end
    S{ind1}=pts(member_tot,:); % 'member_tot' contains the indices for all the connected pixels, here create a separate list for them. 
    pts(member_tot,:)=[]; % Delete elements that have already been taken into account. This speeds up the sorting.
    Mconnect(member_tot,:)=[];
    Mconnect(:,map2rowvec(member_tot))=[];
    ind1=ind1+1;
end
end

function coords = calculate_coordinates(S,minlistlength)

hits=0;
coords=[];
for ind1=1:length(S)
    if length(S{ind1})>minlistlength % checks whether each list contains a number of pixels that is reasonable large
                                     % this helps filter out single pixel "defects" with large illumination value on them, which may result in false hits
        hits=hits+1;
        coords(hits,:)=[0 0];
        for ind2=1:length(S{ind1})
            coords(hits,:)=coords(hits,:)+S{ind1}(ind2,3)*S{ind1}(ind2,1:2);
        end
        coords(hits,:)=coords(hits,:)/sum(S{ind1}(:,3));
    end
end
end
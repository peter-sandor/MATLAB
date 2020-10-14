function varargout = centroid_FF(varargin)
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
else
    pic=varargin{1};
    thrs=varargin{2}.thrs;
%     radius=varargin{2}.radius; 
    minlistlength=varargin{2}.minlistlength;
%     filtersize=varargin{2}.filtersize;
%     max_coord_diff=varargin{2}.max_coord_diff;
end

N=size(pic);
ind1=1;

ind3=1;
offset=0;
while ind1<=N(1)
    ind2=1;
    while ind2<=N(2)
        list_to_check=[];
        list_of_found=[];
        if pic(ind1,ind2)>thrs
            [pic,list_of_found_new,list_to_check_new]=check_neighbors([ind1 ind2 pic(ind1,ind2)],thrs,pic,N,list_of_found,list_to_check);
            pic(ind1,ind2)=0;
            list_of_found=list_of_found_new;
            list_to_check=list_to_check_new;
            while ~isempty(list_to_check)
                [pic,list_of_found_new,list_to_check_new]=check_neighbors(list_to_check(1,:),thrs,pic,N,list_of_found,list_to_check);
                list_of_found=list_of_found_new;
                list_to_check=consolidate_list(list_to_check_new);
            end
            S{ind3}=list_of_found;
            offset=max(list_of_found(:,2))-min(list_of_found(:,2));
            ind3=ind3+1;
        else
            offset=1;
        end
        ind2=ind2+offset;
    end
    ind1=ind1+1;
end
if exist('S')==1
    coords=calculate_coordinates(S,minlistlength);
else
    coords=[];
end
varargout{1}=coords; % return with a list of coordinates
end

function [pic_out list_of_found_new list_to_check_new] = check_neighbors(pixel,thrs,pic_in,N,list_of_found,list_to_check)
% N=size(pic);
pic_out=pic_in;
list_of_found_new=[list_of_found; pixel];
list_to_check_new=list_to_check;
if size(list_to_check_new,1)>1
    list_to_check_new(1,:)=[];
else
    list_to_check_new=[];
end

if pixel(1)==1
    range1=[0 1];
elseif pixel(1)==N(1)
    range1=[0];
else
    range1=[0 1];
end

if pixel(2)==1
    range2=[0 1];
elseif pixel(2)==N(2)
    range2=[-1 0];
else
    range2=[-1 0 1];
end

for ind1=1:length(range1)
    for ind2=1:length(range2)
        if pic_out(pixel(1)+range1(ind1),pixel(2)+range2(ind2))>thrs && (range1(ind1)~=0 || range2(ind2)~=0)
            list_to_check_new=cat(1,list_to_check_new,[pixel(1)+range1(ind1) pixel(2)+range2(ind2) pic_out(pixel(1)+range1(ind1),pixel(2)+range2(ind2))]);
            pic_out(pixel(1)+range1(ind1),pixel(2)+range2(ind2))=0;
        end
    end
end
end

function list_out = consolidate_list(list_in)
ind1=1;
% N=size(list_in,1);
list_out=list_in;
while ind1<size(list_out,1)
    indices1=findvec(list_out(1:ind1-1,:),list_out(ind1,:));
    indices2=findvec(list_out(ind1+1:end,:),list_out(ind1,:));
    indices=logical([indices1; 0; indices2]);
    list_out(indices,:)=[];
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
        for ind2=1:size(S{ind1},1)
            coords(hits,:)=coords(hits,:)+S{ind1}(ind2,3)*S{ind1}(ind2,1:2);
        end
        coords(hits,:)=coords(hits,:)/sum(S{ind1}(:,3));
    end
end
end
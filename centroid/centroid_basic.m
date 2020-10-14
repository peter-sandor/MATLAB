function varargout = centroid_basic(varargin)
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

if nargin==0
    [FileName,PathName] = uigetfile({'*.jpg','JPEG image';'*.png','PNG image';'*.bmp','BMP image';'*.*','All Files (*.*)'},'Pick an image file');
    path_file=strcat(PathName,FileName);
    pic=monopic(imread(path_file));
    thrs=10;
    radius=2; % distance based on which it is determined whether two pixels are topologically connected or not (are they close enough to each other)
    minlistlength=4; % this number sets the threshold for a check later on in the code that tells whether each list contains a number of pixels that is reasonable large
    filtersize=60;
else
    pic=varargin{1};
    thrs=varargin{2}.thrs;
    radius=varargin{2}.radius; % distance based on which it is determined whether two pixels are topologically connected or not (are they close enough to each other)
    minlistlength=varargin{2}.minlistlength; % this number sets the threshold for a check later on in the code that tells whether each list contains a number of pixels that is reasonable large
    filtersize=varargin{2}.filtersize;
end

pts=thrs_img(pic,thrs,radius,minlistlength,filtersize); % threshold image and enumerate the points surviving the process --> 'pts'
if size(pts,1)>1    % if the list is not empty, we start sorting its elements into lists, which will corresponds to the different hits on the detector
    S=cluster_pixels(pts,radius,minlistlength);
    % Done with creating the lists, next step is to see if they are 'real'
    % hits. In case yes, their center coordinates are calculated.
    coords = calculate_coordinates(S,minlistlength);
else
    coords=[];
end
varargout{1}=coords; % return with a list of coordinates
end
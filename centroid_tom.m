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
tic
for k=1:100
%pic(pic<=thrs)=0; % thresholding is done here
[y,x]=find(pic>=thrs);
tom1=circshift([y,x],1);
tom=sum([tom1(:,1)-y,tom1(:,2)-x].^2,2);
index=(tom<=10);
%fft2(pic);
end
tom
size(tom)
size([x,y])
[x(index),y(index)]

toc
% if max(max(pic))>0 % if there are pixels with nonzero values still, we start sorting them into the lists
%     pic=double(pic);
%     [y,x]=find(pic);
%     pts(:,1)=x;
%     pts(:,2)=y;
%     M=size(pts,1);
%     % Now create connection matrix with size M x M, logical type. Two pixels are connected
%     % if their distance is smaller than the value of 'radius'.
%     for ind1=1:M
%         Mconnect(:,ind1)=logical(((pts(ind1,1)-pts(:,1)).^2+(pts(ind1,2)-pts(:,2)).^2)<=radius^2);
%     end
%     
%     % The next part uses the connection matrix to sort the pixels.
%     ind2=1;
%     while ~isempty(Mconnect)
%     % find first-order connections
%     member0=logical(Mconnect(:,1)); % index vector (column)
%     % find second-order connections
%     member1=any(Mconnect(member0,:),1);
%     member1=map2colvec(member1); % index vector (column)
%     member_diff=xor(member0,member1);
%     Mconnect(member0,:)=false;
% 	Mconnect(:,map2rowvec(member0))=false;
%     member_tot=member1;
%     % check and follow up on higher-order connections
%     while any(member_diff)
%         member0=member1;
%         member1=any(Mconnect(member_diff,:),1);
%         member1=map2colvec(logical(member1)); % index vector (column)
%         member_diff=xor(member0,member1);
%         member_tot=or(member_tot,member1);
%         Mconnect(member0,:)=false;
%         Mconnect(:,map2rowvec(member0))=false;
%     end
%     S{ind2}=pts(member_tot,:); % 'member_tot' contains the indices for all the connected pixels, here create a separate list for them. 
%     pts(member_tot,:)=[];
%     Mconnect(member_tot,:)=[];
%     Mconnect(:,map2rowvec(member_tot))=[];
%     ind2=ind2+1;
%     end
%     
%     % Done with creating the lists, next step is to see if they are 'real'
%     % hits. In case yes, their center coordinates are calculated.
%     hits=0;
%     cntr=[];
%     for ind1=1:length(S)
%         if length(S{ind1})>minlistlength % checks whether each list contains a number of pixels that is reasonable large
%                                          % this helps filter out single pixel "defects" with large illumination value on them, which may result in false hits
%             hits=hits+1;
%             cntr(hits,:)=[0 0];
%             weights=zeros([length(S{ind1}) 1]);
%             for ind2=1:length(S{ind1})
%                 weights(ind2)=pic(S{ind1}(ind2,2),S{ind1}(ind2,1));
%                 cntr(hits,:)=cntr(hits,:)+weights(ind2)*S{ind1}(ind2,:);
%             end
%             cntr(hits,:)=cntr(hits,:)/sum(weights);
%         end
%     end
% else cntr=[];
%     hits=0;
% end
% varargout{1}=cntr; % return with a list of coordinates
% end
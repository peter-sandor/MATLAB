function h = imagescP(varargin)
%IMAGESC Scale data and display as image.
%   IMAGESC(...) is the same as IMAGE(...) except the data is scaled
%   to use the full colormap.
%   
%   IMAGESC(...,CLIM) where CLIM = [CLOW CHIGH] can specify the
%   scaling.
%
%   See also IMAGE, COLORBAR, IMREAD, IMWRITE.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.11.4.5 $

clim = [];
switch (nargin),
  case 0,
    hh = image('CDataMapping','scaled');
  case 1,
    hh = image(varargin{1},'CDataMapping','scaled');
  case 3,
    hh = image(varargin{:},'CDataMapping','scaled');
  otherwise,

    % Determine if last input is clim
    if isequal(size(varargin{end}),[1 2])
      str = false(length(varargin),1);
      for n=1:length(varargin)
        str(n) = ischar(varargin{n});
      end
      str = find(str);
      if isempty(str) || (rem(length(varargin)-min(str),2)==0),
        clim = varargin{end};
        varargin(end) = []; % Remove last cell
      else
        clim = [];
      end
    else
      clim = [];
    end
    hh = image(varargin{:},'CDataMapping','scaled');
end

% Get the parent Axes of the image
cax = ancestor(hh,'axes');
set(cax,'ydir','normal');
if ~isempty(clim),
  set(cax,'CLim',clim)
elseif ~ishold(cax),
  set(cax,'CLimMode','auto')
end

if nargout > 0
    h = hh;
end

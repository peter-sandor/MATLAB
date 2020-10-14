function out = autoconv(in,varargin)
if nargin == 1 % nargin checks for the total number of arguments
    out = conv(in,in);
elseif nargin == 2
    out = conv(in,in,varargin{1});
else
    disp('Wrong number of arguments (have to be either 1 or 2)!');
end
end
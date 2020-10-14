function out = MCPmask(in,center,varargin)
% MCP sensitivity mask based on imaging light ions of alpha-terpinene on
% 02/19/14
% for the construction of the mask explicitly, see
% D:\users\peter\matlab_processing\create_sens_mask2.m
if nargin==2
    load 'D:\users\peter\matlab_processing\mask_2014feb.mat';
    mask=mask_2014feb;
elseif nargin==3;
    mask=varargin{1};
end
if findvec(size(in),size(mask))
    out = in.*circshift(mask,flipdim(center-[190 200],2)); % center: [x,y]=[column,row]
else
    disp('dimensions don''t match');
    out=[];
end
function struct2var(S,varargin)
% This function takes struct 'S' and creates variables using its field names and field values
% Optional argument: workspace name in which the variables are assigned to.

if nargin==1
    wspace = 'base';
else
    wspace = varargin{1};
end
names = fieldnames(S);
for ind1=1:length(names)
    assignin(wspace,names{ind1},getfield(S,names{ind1}));
end
end
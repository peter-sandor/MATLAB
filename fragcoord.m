function [Frags, Coords] = FragCoord(filename)

% this code takes the results of a coincidence measurement ('coords.txt')
% and writes the coordinates of hits corresponding to the different
% fragments into separate files ('fragment_N.txt' where 'N' is the ID for a specific fragment).

% Coods=dlmread(filename);
Coods=load(filename);
[a,b,c]=fileparts(filename);
if ~isempty(a)
    a=[a '\'];
end
Frags=unique(Coods(:,1));
Frags(Frags==0)=[];
for ind1=1:length(Frags)
    dlmwrite([a 'fragment_' num2str(Frags(ind1)) '.txt'],Coods(Coods(:,1)==Frags(ind1),2:3),'\t');
end
end
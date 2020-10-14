function clr(filename)

s=whos('-file',filename);
for ind1=1:length(s)
    toclear=s(ind1).name;
    evalin('base',['clear(''' toclear ''')'])
%     evalin('base',['clear(' toclear ')'])
%     evalin('base',['clear(s(' num2str(ind1) ').name)']);
end
end
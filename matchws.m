function varargout = matchws(filename)
c=evalin('base','who');
s=whos('-file',filename);
flag=0;
ind3=1;
for ind1=1:length(c)
    for ind2=1:length(s)
        if strcmp(c{ind1},s(ind2).name)
            varargout{1}{ind3}=c{ind1};
            disp(c{ind1});
            ind3=ind3+1;
            flag=1;
        end
    end
end
if flag==0
    varargout{1}=[];
    disp('No match has been found.');
end
end

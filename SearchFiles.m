function ind_out = SearchFiles(files,Name)

N=length(files);
ind_out=[];
for ind1=1:N
    if strcmp(files(ind1).name,Name)
        ind_out=[ind_out ind1];
    end
end

end
    
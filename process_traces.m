files=dir('trace*');
fid=fopen('trace_counts.txt','w');
N=length(files);
counter=0;
for ind1=1:N
    hits(ind1)=tracefilter(files(ind1).name,133);
    if hits(ind1)>0
        fprintf(fid,'%s \n',[files(ind1).name '    ' num2str(hits(ind1))]);
        counter=counter+1;
    end
    disp([num2str(ind1) '/' num2str(N)]);
end
fclose(fid);
disp(['traces with nonzero hits:' num2str(counter) '/' num2str(N)]);
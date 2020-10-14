files=dir('trace*');
fid_tr=fopen('trc_counts.txt','w');
fid_im=fopen('img_counts.txt','w');
fid_both=fopen('both_counts.txt','w');
N=length(files);
counter_both=0;
counter_im=0;
counter_tr=0;
for ind1=1:N
    hit_trace(ind1)=tracefilter(files(ind1).name,122,[2500 2600],1);
    hit_image(ind1)=imfilter2(['pic' files(ind1).name(6:end) '.bin'],40,[170 200;180 210]);
    if hit_image(ind1)>0
        if hit_trace(ind1)>0
            fprintf(fid_both,'%s \n',[files(ind1).name(6:end) '    ' num2str(1)]);
            counter_both=counter_both+1;
        else fprintf(fid_im,'%s \n',[files(ind1).name(6:end) '    ' num2str(1)]);
            counter_im=counter_im+1;
        end
    elseif hit_trace(ind1)>0
        fprintf(fid_tr,'%s \n',[files(ind1).name(6:end) '    ' num2str(1)]);
        counter_tr=counter_tr+1;
    end
    if mod(ind1,100)==0
        disp([num2str(ind1) '/' num2str(N)]);
    end
end
fclose(fid_tr);
fclose(fid_im);
fclose(fid_both);
disp(['both:' num2str(counter_both) '/' num2str(N)]);
disp(['trace only:' num2str(counter_tr) '/' num2str(N)]);
disp(['image only:' num2str(counter_im) '/' num2str(N)]);
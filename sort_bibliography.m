function varargout = sort_bibliography(texname,bibname)
% This code parses a LaTex file (.tex) and collects the citation handles.
% Then it takes the entries corresponding to the handles from a separate
% file and writes them to a third file in order of appearance.
str_cttn = get_cite_order(texname);
ind_commas=strfind(str_cttn,',');
if ~isempty(ind_commas)
    N_commas=length(ind_commas);
    citation(1).ID=str_cttn(1:ind_commas(1)-1);
    ind_cit=2;
    for ind1=2:N_commas
        % make sure the new new entry is not present in the list compiled
        % so far
        str_new=str_cttn(ind_commas(ind1-1)+1:ind_commas(ind1)-1);
        flag_new=1;
        for ind3=1:length(citation)
            if strcmp(str_new,citation(ind3).ID)
                flag_new=0;
            end
        end
        if flag_new
            citation(ind_cit).ID=str_new;
            ind_cit=ind_cit+1;
        end
    end
    str_new=str_cttn(ind_commas(end)+1:end);
    flag_new=1;
    for ind3=1:length(citation)
        if strcmp(str_new,citation(ind3).ID)
            flag_new=0;
        end
    end
    if flag_new
        citation(ind_cit).ID=str_new;
    end
end
Nid=length(citation);
fid1=fopen(bibname,'r');
fid2=fopen([bibname(1:end-4) '_sorted' bibname(end-3:end)],'w');

for ind_ID = 1:Nid
    frewind(fid1);
    content=[];
    ind_C=1;
    lastpos=ftell(fid1);
    % find a specific citation and copy its contents
    temp=fgetl(fid1);
%         ind_bib1=strfind(temp,'\bibitem')+9;
%         ind_bib2=strfind(temp(ind_bib1:end),'}')+ind_bib1-2;

    while (isempty(strfind(temp,'\bibitem'))) || (~strcmp(temp(strfind(temp,'\bibitem')+9:strfind(temp(strfind(temp,'\bibitem')+9:end),'}')+strfind(temp,'\bibitem')+9-2),citation(ind_ID).ID))
        lastpos=ftell(fid1);
        if feof(fid1)
            flag1=1;
            break;
        else flag1=0;
        end
        temp=fgetl(fid1);
    end
    if flag1==0
        fseek(fid1,lastpos,'bof');
        citation(ind_ID).content{ind_C}=fgetl(fid1);
        temp=fgetl(fid1);
        if ind_ID==80
            1;
        end
        while isempty(strfind(temp,'\bibitem')) && ~feof(fid1)
            ind_C=ind_C+1;
            citation(ind_ID).content{ind_C}=temp;
            temp=fgetl(fid1);
        end
        Nc=ind_C-1;
        for ind2=1:Nc
            fprintf(fid2,'%s\n',citation(ind_ID).content{ind2});
        end
        fprintf(fid2,'\n');
    end
end
varargout{1}=citation;
fclose(fid1);
fclose(fid2);
end
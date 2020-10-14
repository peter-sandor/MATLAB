function date_out = folder2date(foldername)
ind_bslash=strfind(foldername,'\');
if length(ind_bslash)>=2;
    if ind_bslash(end)==length(foldername)
        date_out=StringReplace(foldername(ind_bslash(end-1)+1:ind_bslash(end-1)+10),'_','-');
    else
        date_out=StringReplace(foldername(ind_bslash(end)+1:ind_bslash(end)+10),'_','-');
    end
else
    date_out=[];
    disp('Can not convert foldername to date.')
end
end
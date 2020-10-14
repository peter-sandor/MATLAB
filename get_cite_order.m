function str_cttn = get_cite_order(fname)
fid=fopen(fname,'r');
if fid~=-1
    str_cttn = [];
    while ~feof(fid)
       str1=fgetl(fid);
       index1=strfind(str1,'\cite{');
       Nind=length(index1);
       if ~isempty(index1)
           index1=index1+6;
           if Nind>=2
               for ind1=1:Nind-1
                temp=strfind(str1(index1(ind1):index1(ind1+1)-1),'}');
                index2(ind1)=temp(1)+index1(ind1)-2;
                str_cttn=[str_cttn ',' str1(index1(ind1):index2(ind1))];
               end
           end
           ind1=Nind;
           temp=strfind(str1(index1(ind1):length(str1)),'}');
           index2(ind1)=temp(1)+index1(ind1)-2;
           str_cttn=[str_cttn ',' str1(index1(ind1):index2(ind1))];
       end
    end
    fclose(fid);
end
str_cttn=str_cttn(2:end);
end
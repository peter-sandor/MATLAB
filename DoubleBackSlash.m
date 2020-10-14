function strout = DoubleBackSlash(strin)
index=strfind(strin,'\');
N=length(index);
strout=[strin(1:index(1)) '\'];
for ind1=1:N-1
    strout=[strout strin((index(ind1)+1):index(ind1+1)) '\'];
end
strout=[strout strin(index(N)+1:length(strin))];
end
    
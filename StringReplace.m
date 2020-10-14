function str_out = StringReplace(str_in,lookfor,replacewith)

Nlook=length(lookfor);
Nrepl=length(replacewith);
str_out=str_in;
ind_str=strfind(str_out,lookfor);
while ~isempty(ind_str)
    str_out=[str_out(1:ind_str(1)-1) replacewith str_out(ind_str(1)+Nlook:end)];
    ind_str=strfind(str_out(ind_str(1)+Nrepl:end),lookfor)+ind_str(1)+Nrepl-1;
end
end
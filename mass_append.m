function mass_append(string_to_append)
c=evalin('base','who');
N=length(c);
for ind1=1:N
    evalin('base',[c{ind1} string_to_append '=' c{ind1} ';']);
    evalin('base',['clear(''' c{ind1} ''');']);
    clear(c{ind1});
end
end
function string_out = put_array_ops(string_in)

string_out=string_in;
operators=['*/^'];

for ind2=1:length(operators)
    inds1=strfind(string_out,operators(ind2));
    for ind1=length(inds1):-1:1
        if string_out(inds1(ind1)-1)=='.'
            inds1(ind1)=[];
        end
    end
    while ~isempty(inds1)
        N=length(string_out);
        string_out(inds1(end)+1:N+1)=string_out(inds1(end):N);
        string_out(inds1(end):inds1(end)+1)=['.' operators(ind2)];
        inds1=strfind(string_out,operators(ind2));
        for ind1=length(inds1):-1:1
            if string_out(inds1(ind1)-1)=='.'
                inds1(ind1)=[];
            end
        end
    end
end

end
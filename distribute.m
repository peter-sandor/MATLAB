function out = distribute(Nballs,Nslots)
out=[];
vec=Nballs;
oldrows=[vec zeros([1 Nslots-length(vec)])];
out=[out; oldrows];
while size(oldrows,1)>=1
    newrows=[];
    for ind1=1:size(oldrows,1)
        vec=oldrows(ind1,oldrows(ind1,:)~=0);
        if (length(vec)<=Nslots) && isempty(findvec(vec,ones([1,Nballs])))
            getindex=(1:length(vec)).*(vec>1);
            getindex=getindex(getindex>0);
            for ind3=1:length(getindex)
                putindex=unique(((getindex(ind3)+1):length(vec)).*(vec((getindex(ind3)+1):min(end,Nslots))<(vec(getindex(ind3))-1)));
                if max(putindex)==0;
                    putindex=[];
                end
                if length(vec)<Nslots
                    putindex=unique([putindex 0]);
                    Nput=length(putindex);
                    for ind2=1:Nput
                        newrows=[newrows; makestep(oldrows(ind1,:),getindex(ind3),putindex(ind2))];
                    end
                elseif length(vec)==Nslots
                    putindex(putindex==0)=[];
                    Nput=length(putindex);
                    for ind2=1:Nput
                        newrows=[newrows; makestep(oldrows(ind1,:),getindex(ind3),putindex(ind2))];
                    end
                end
            end
        end
    end
    out=[out; newrows];
    oldrows=uniquevec(newrows);
end
out=uniquevec(sortrows(out));

	function rowout = makestep(rowin,getind,putind)
        vec=rowin(rowin~=0);
        vec(getind)=vec(getind)-1;
        if isempty(putind) || putind==0 || (putind==length(vec)+1)
            vec=[vec 1];
        else
            vec(putind)=vec(putind)+1;
        end
        rowout=sort([vec zeros([1 length(rowin)-length(vec)])],'descend');
    end

end
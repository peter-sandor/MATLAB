function out=xcorrP(in1,in2)
N=max(length(in1),length(in2));
M=min(length(in1),length(in2));
in1=map2colvec(in1);
in2=map2colvec(in2);
out=zeros([N+M-1 1]);
norm=sum(in1)*sum(in2);
if length(in2)>length(in1)
    temp=in1;
    in1=in2;
    in2=temp;
    flag=1;
    for ind1=1:N+M-1
        out(ind1)=sum(conj(in1(max(N-ind1+1,1):min(N,N-ind1+M))).*in2(max(ind1-N+1,1):min(ind1,M)));%/norm;
%   disp(ind1)
    end
    out=flipdim(out,1);
else    for ind1=1:N+M-1
            out(ind1)=sum(in1(max(N-ind1+1,1):min(N,N-ind1+M)).*conj(in2(max(ind1-N+1,1):min(ind1,M))))/norm;
        end
%     disp(ind1)
end
end
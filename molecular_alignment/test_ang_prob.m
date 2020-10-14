Ylm=simresult.Ylm_red;
amps=zeros([size(Ylm,2) size(Ylm,3) size(simresult.amps0,4)]);
for indJstate=1:size(simresult.amps0,4)-1
    for indJ=1:size(Ylm,2)
        for indM=1:min([indJ size(Ylm,3)])
            if indJ==indJstate
                if indM==1
                    amps(indJ,indM,indJstate)=1;
                else
                    amps(indJ,indM,indJstate)=2;
                end
            end
        end
    end
end
Ylm_sum=zeros([size(Ylm,1) size(simresult.amps0,4)]);
for ind6=1:size(simresult.amps0,4)
    Ylm_sum(:,ind6) = simresult.JM_probs(ind6)*map2colvec(sum((sum(Ylm.*permute(extend(squeeze(amps(:,:,ind6)),size(Ylm,1)),[3 1 2]),2)).^1,3));
end
Ylm_sum1=squeeze(sum(Ylm_sum,2));
Ylm_norm=Ylm_sum1./sum(Ylm_sum1.*map2colvec(sin(simresult.theta)));
figure;imagescP(Ylm_sum)
figure;plot(Ylm_sum)
figure;plot(Ylm_norm)
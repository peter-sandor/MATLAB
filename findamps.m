function ampsout = findamps(ints,allamps)
eps=1e-12;
ints=map2colvec(ints);
allamps=map2colvec(allamps);
temp=amp2int(allamps);
ampsout=zeros(size(ints));
for ind1=1:length(ints)
    for ind2=1:length(allamps)
        if abs(ints(ind1)-temp(ind2))<=eps;
            ampsout(ind1)=allamps(ind2);
            break;
        end
    end
end
end

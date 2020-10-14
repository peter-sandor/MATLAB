function P = ProbIonDet2(lambda,electronstats,ionstats)
% This function calculates the probability of generating ion+electron pairs in the laser field and detecting those with parameters given in
% the matrices 'electronstats' and 'ionstats'. Electronstats only has one row, while each row of 'ionstats' contain values for a given
% fragment type. size(ionstats,1) = number of different fragments types
% lambda is the expectation value of the number of pairs generated,
% assuming Poissonian statistics. (Thesis notation is used.)

% electrontats: first column: occurrence number (k)
%           second column: number of detections (m)
%           third column: detection efficiency  (p)

% ionstats: first column: occurrence numbers    (nj)
%           second column: number of detections (lj)
%           third column: detection efficiences (qj)
%           fourth column: branching ratios     (bj)

P = nchoosek(electronstats(1),electronstats(2)).*electronstats(3).^electronstats(2).*(1-electronstats(3)).^(electronstats(1)-electronstats(2)).*lambda.^electronstats(1).*exp(-lambda);
indices=unique((map2colvec(~(ionstats(:,4)==0)).*map2colvec(1:size(ionstats,1))));
indices(indices==0)=[];
for ind2=1:length(indices)
    P=P.*ionstats(indices(ind2),4)^ionstats(indices(ind2),1)/factorial(ionstats(indices(ind2),1)).*nchoosek(ionstats(indices(ind2),1),ionstats(indices(ind2),2)).*ionstats(indices(ind2),3).^ionstats(indices(ind2),2).*(1-ionstats(indices(ind2),3)).^(ionstats(indices(ind2),1)-ionstats(indices(ind2),2));
end
end
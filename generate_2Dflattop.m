center=[527;399]; % x,y
radius=350;
index_fig=zeros([768 1024]);
for ind1=1:size(index_fig,2)
	for ind2=1:size(index_fig,1)
        if ((ind1-center(1))^2+(ind2-center(2))^2)<=radius^2
            index_fig(ind2,ind1)=1;
        end
	end
end
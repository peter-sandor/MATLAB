function varargout = TOFpeakfind(traces,index_search,index_limits,index_bkg,polarity)

ind_check=1;
if polarity>0
    temp=vec2ind(traces(index_search(1):index_search(2),ind_check)==max(traces(index_search(1):index_search(2),ind_check)));
else
    temp=vec2ind(traces(index_search(1):index_search(2),ind_check)==min(traces(index_search(1):index_search(2),ind_check)));
end
ind_peak=index_search(1)-1+temp(1);
for ind1=1:size(traces,2)
    boxcar(ind1)=map2colvec(squeeze(sum(traces(ind_peak+index_limits(1):ind_peak+index_limits(2),ind1)-mean(traces(index_bkg(1):index_bkg(2),ind1),1))));
end
boxcar=map2colvec(boxcar);
amp=map2colvec(squeeze(traces(ind_peak,:)-traces(ind_peak+index_limits(1),:)));

varargout{1}=ind_peak;
varargout{2}=boxcar;
varargout{3}=amp;%Amplitude?
end
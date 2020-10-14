function out = avgtrace(namestub,pcol,pvalue)
% This code averages traces in the current directory corresponding to a
% specific value in the parameterspace file. This parameter value is
% supplied in the argument: 'pvalue', while 'pcol' is the column in the
% parameterspace file in which the code looks for 'pvalue'.
% 'namestub' is the common beginning of the filenames (e.g.: 'pic_1001.jpg' --> namestub='pic_';)
% Thresholding is applied: values <= thrs are set to zero.
% Typical application: intensity scan: a0+a1*0 --> average traces of
% specific AOM amplitude
% files = dir(strcat('*.',ext));
ext='dat';
Npics = load('NumberOfTraces.txt');
pspace0 = load('paramtersspace.txt');
% index_vector=logical(pspace0(:,pcol)==pvalue);
index_vector=logical(findvec(pspace0(:,[1:pcol-1 pcol+1:end]),pvalue));
bkg_vector=logical(Npics<0);
index_vector=index_vector-bkg_vector;
index_vector(index_vector<0)=0;
filenums=unique((0:length(index_vector)-1).*map2rowvec(index_vector));
bkgnums=unique((0:length(bkg_vector)-1).*map2rowvec(bkg_vector));
filenums(filenums==0)=[];
bkgnums(bkgnums==0)=[];
Nvec=length(filenums);
Nbkg=length(bkgnums);
sumtrace=zeros(size(load([namestub num2str(filenums(1)) '.' ext])));
for ind = 1:Nvec
    curname=[namestub num2str(filenums(ind)) '.' ext];
    if exist(curname)==2
        trace=load(curname);
        trace=trace-trace(1);
    else
        continue;
    end        
    sumtrace = sumtrace + trace;
end;
out=sumtrace/Nvec;
end
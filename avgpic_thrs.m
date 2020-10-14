function out = avgpic_thrs(namestub,ext,pcol,pvalue,thrs)
% This code averages pictures in the current directory corresponding to a
% specific value in the parameterspace file. This parameter value is
% supplied in the argument: 'pvalue', while 'pcol' is the column in the
% parameterspace file in which the code looks for 'pvalue'.
% 'namestub' is the common beginning of the filenames (e.g.: 'pic_1001.jpg' --> namestub='pic_';)
% Thresholding is applied: values <= thrs are set to zero.
% Typical application: intensity scan: a0+a1*0 --> average images of
% specific AOM amplitude
% The code handles binary images (.bin) and JPEGs as well (.jpg) --> 'ext'
% load detD;
% files = dir(strcat('*.',ext));
Npics = load('NumberOfPictures.txt');
pspace0 = load('paramtersspace.txt');
% index_vector=logical(pspace0(:,pcol)==pvalue);
index_vector=logical(findvec(pspace0(:,[1:pcol-1 pcol+1:end]),pvalue));
bkg_vector=logical(Npics<0);
index_vector=index_vector-bkg_vector;
index_vector(index_vector<0)=0;
filenums=unique((1:length(index_vector)).*map2rowvec(index_vector));
bkgnums=unique((1:length(bkg_vector)).*map2rowvec(bkg_vector));
filenums(filenums==0)=[];
bkgnums(bkgnums==0)=[];
Nvec=length(filenums);
Nbkg=length(bkgnums);
if strcmp(ext,'jpg')==1
	sumpic=zeros(size(monopic(imread([namestub num2str(filenums(1)) '.' ext]))));
elseif strcmp(ext,'bin')==1
	sumpic=zeros(size(loadbinaryimage([namestub num2str(filenums(1)) '.' ext])));
end

for ind = 1:Nvec
    curname=[namestub num2str(filenums(ind)) '.' ext];
    if exist(curname)==2
        if strcmp(ext,'jpg')==1
            pic=monopic(imread(curname));
        elseif strcmp(ext,'bin')==1
            pic=loadbinaryimage(curname);
        end
        pic(pic<=thrs)=0; % threshold image
    else
        continue;
    end        
    sumpic = sumpic + pic;
end;
out=sumpic/Nvec;
end



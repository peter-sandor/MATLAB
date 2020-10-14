function allcoords = imgreconst(namestub,ext,thrs,pcol,pvalue)
% This code get the coordinates of hits on the images in the current directory corresponding to a
% specific value in the parameterspace file. This parameter value is
% supplied in the argument: 'pvalue', while 'pcol' is the column in the
% parameterspace file in which the code looks for 'pvalue'.
% 'namestub' is the common beginning of the filenames (e.g.: 'pic_1001.jpg' --> namestub='pic_';)
% Typical application: intensity scan: a0+a1*0 --> average images of
% specific AOM amplitude
% The code handles binary images (.bin) and JPEGs as well (.jpg) --> 'ext'
% load detD;
% files = dir(strcat('*.',ext));
pspace0 = load('paramtersspace.txt');
Npics = load('NumberOfPictures.txt');
bkg_vector=logical(Npics<0);
% index_vector=logical(pspace0(:,pcol)==pvalue);
index_vector=logical(findvec(pspace0(:,[1:pcol-1 pcol+1:end]),pvalue));
index_vector=index_vector-bkg_vector;
index_vector(index_vector<0)=0;
filenums=unique((1:length(index_vector)).*map2rowvec(index_vector));
filenums(filenums==0)=[];
Nvec=length(filenums);
allcoords=[];
for ind = 1:Nvec
    curname=[namestub num2str(filenums(ind)) '.' ext];
    if exist(curname)==2
        if strcmp(ext,'jpg')==1
            pic=monopic(imread(curname));
        elseif strcmp(ext,'bin')==1
            pic=loadbinaryimage(curname);
        end
    else
        continue;
    end      
    [a,b]=getshots(pic,thrs);
    disp(['image:' num2str(ind) '/' num2str(Nvec)]);
    allcoords = cat(1,allcoords,b);
end
end
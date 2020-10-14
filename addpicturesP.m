function out = addpicturesP(pcol,pvalue)
% extension si file extension
% lastmain is last number of mainpictures
load detD;
files = dir('image*.jpg');
Npics = load('NumberOfPictures.txt');
pspace0 = load('paramtersspace.txt');
pspace=pspace0;
pspace(Npics<0,:)=[]; % delete rows from parameterspace array which correspond to background image
index_vector=logical(pspace0(:,pcol)==pvalue);
vec1=1 : length(files);
sumpic=zeros(size(monopic(imread(files(1).name))));
for ind = vec1(index_vector)
    pic = monopic(imread(files(ind).name));
    sumpic = sumpic + pic;
end;
pic=sumpic/length(vec1(index_vector));
bkgpic=zeros(size(sumpic));
for ind=vec1(logical(Npics<0))
    pic=monopic(imread(files(ind).name));
    bkgpic=bkgpic+pic;
end
bkgpic=bkgpic/length(vec1(logical(Npics<0)));
save sumpic.mat sumpic;
save bkgpic.mat bkgpic
out = (sumpic - bkgpic)./detD;
end



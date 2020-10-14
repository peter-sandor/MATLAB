function out = avgpic(ext,pcol,pvalue)
% This code averages pictures in the current directory corresponding to a
% specific value in the parameterspace file. This parameter value is
% supplied in the argument: 'pvalue', while 'pcol' is the column in the
% parameterspace file in which the code looks for 'pvalue'.
% Typical application: intensity scan: a0+a1*0 --> average images of
% specific AOM amplitude
% The code handles binary images (.bin) and JPEGs as well (.jpg) --> 'ext'
% load detE;
files = dir(strcat('*.',ext));
Npics = load('NumberOfTraces.txt');
pspace0 = load('paramtersspace.txt');
pspace=pspace0;
pspace(Npics<0,:)=[]; % delete rows from parameterspace array which correspond to background image
index_vector=logical(pspace(:,pcol)==pvalue);
vec1=1 : length(files);
if strcmp(ext,'jpg')==1
	sumpic=zeros(size(monopic(imread(files(1).name))));
elseif strcmp(ext,'bin')==1
	sumpic=zeros(size(loadbinaryimage(files(1).name)));
end
for ind = vec1(index_vector)
    if strcmp(ext,'jpg')==1
        pic = monopic(imread(files(ind).name));
    elseif strcmp(ext,'bin')==1
        pic = loadbinaryimage(files(ind).name);
    end
    sumpic = sumpic + pic;
end;
pic=sumpic/length(vec1(index_vector));
bkgpic=zeros(size(sumpic));
for ind=vec1(logical(Npics<0))
    if strcmp(ext,'jpg')==1
        pic=monopic(imread(files(ind).name));
    elseif strcmp(ext,'bin')==1
        pic = loadbinaryimage(files(ind).name);
    end
    bkgpic=bkgpic+pic;
end
bkgpic=bkgpic/length(vec1(logical(Npics<0)));
save sumpic.mat sumpic;
save bkgpic.mat bkgpic
% out = (sumpic - bkgpic).*detE;
out = (sumpic - bkgpic);
end



function [out,varargout] = avgpic(ext,pvec)
% This code averages pictures in the current directory corresponding to a
% specific vector of values in the parameterspace file. This vector is
% supplied in the argument: 'pvec'. The goal is to average over different
% values of dummy variables. The code assumes that the paramater of the
% largest index (in case of {a0,a1,a2}, it's a2) is the dummy variable -->
% last column of 'paramtersspace' file.
% Typical application: combined intensity and pump-probe scan --> average of images of
% specific AOM amplitude and specific delay.
% The code handles binary images (.bin) and JPEGs as well (.jpg) --> 'ext'
load detD;
files = dir(strcat('*.',ext));
Npics = load('NumberOfPictures.txt');
pspace0 = load('paramtersspace.txt');
pspace=pspace0;
pspace(Npics<0,:)=[]; % delete rows from parameterspace array which correspond to background image
Ndummy=size(pspace,2)-length(pvec)
index_vector=logical(zeros([size(pspace,1) 1]));
for ind1=1:size(pspace,1)
    index_vector(ind1)=logical(isequal(pspace(ind1,[1:(end-Ndummy)]),pvec));
end
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
varargout{1}=bkgpic;
save sumpic.mat sumpic;
save bkgpic.mat bkgpic
out = (sumpic - bkgpic)./detD;
end



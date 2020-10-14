function varargout = copy_tex_figures(varargin)
% This code extracts source files for graphics objects from a latex
% document (.tex), and backs them up to a user-defined path.

% Input arguments: filename -- .tex file to process
%                  copy_to -- target directory; this is an optional
%                               argument, if not supplied, the code just
%                               extracts the paths for graphics

% Output arguments: filepath -- cell array containing the paths of the
%                               graphics files

filename=varargin{1};
fileID=fopen(filename);
str_to_find1='\includegraphics'; % Graphics entries are identified by the command '\includegraphics'.
Nstr1=length(str_to_find1);
filepath=[];

ind1=1;
tline = fgetl(fileID);
index1=strfind(tline,str_to_find1);
if ~isempty(index1);
    tline2=tline(index1+Nstr1:end);
    if tline2(1)=='[';
        index2=strfind(tline2,']');
        index3=min(strfind(tline2,'}'));
        str1=tline2(index2+2:index3-1);
    else
        index3=min(strfind(tline2,'}'));
        str1=tline2(2:index3-1);
    end
    index4=strfind(str1,'/');
    for ind2=1:length(index4)
        str1(index4(ind2))='\';
    end
    filepath{ind1}=str1;
    ind1=ind1+1;
end
while ischar(tline)
    tline = fgetl(fileID);
    index1=strfind(tline,str_to_find1);
    if ~isempty(index1);
        tline2=tline(index1+Nstr1:end);
        if tline2(1)=='[';
            index2=strfind(tline2,']');
            index3=min(strfind(tline2,'}'));
            str1=tline2(index2+2:index3-1);
        else
            index3=min(strfind(tline2,'}'));
            str1=tline2(2:index3-1);
        end
        index4=strfind(str1,'/');
        for ind2=1:length(index4)
            str1(index4(ind2))='\';
        end
        filepath{ind1}=str1;
        ind1=ind1+1;
    end
end
fclose(fileID);

Nfiles=length(filepath);
disp(['Found ' num2str(Nfiles) ' reference(s) to figures in .tex document.']);
ind3=0;
if ~isempty(filepath) && nargin>1
    copy_to=varargin{2};
    if exist(copy_to)~=7
        mkdir(copy_to); % If predefined destination does not exist, create it.
    end
    for ind1=1:Nfiles
        if exist(filepath{ind1})==2
            copyfile(filepath{ind1},copy_to); % Copy graphics file to predefined destination.
            if exist([filepath{ind1}(1:end-3) 'fig'])==2
                copyfile([filepath{ind1}(1:end-3) 'fig'],copy_to); % If exists, copy the corresponding .fig file as well.
            end
            ind3=ind3+1;
        end
    end
end
disp([num2str(ind3) ' files were copied.']);
varargout{1}=filepath;
end
function varargout = secretary(varargin)

if nargin==0
    regdir='D:\users\peter\VMI\';
    regfile='datasets.xlsx';
    basepath='.\';
    filename='comments.txt';
elseif nargin==1
    regdir='D:\users\peter\VMI\';
    regfile='datasets.xlsx';
    basepath=varargin{1};
    filename='comments.txt';
elseif nargin==2
    [regdir,regname,regext]=fileparts(varargin{2});
    regdir=[regdir '\'];
    regfile=[regname regext];
    basepath=varargin{1};
    filename='comments.txt';
elseif nargin==3
    [regdir,regname,regext]=fileparts(varargin{2});
    regdir=[regdir '\'];
    regfile=[regname regext];
    basepath=varargin{1};
    filename=varargin{3};
end
origpath=pwd;

cd(basepath)
files=subdir(filename);
Nfiles=length(files);
disp(['Found ' num2str(Nfiles) ' instance(s) of "' filename '"'])
toxls=[];
maxlength=0;
for ind1=1:Nfiles
	fileID=fopen(files(ind1).name);
	comments=textscan(fileID,'%s');
    fclose(fileID);
    folders{ind1}=files(ind1).name(1:max(strfind(files(ind1).name,filename))-1);
    if length(folders{ind1})>maxlength
        maxlength=length(folders{ind1});
    end
	toxls{ind1,1}=files(ind1).date;
    toxls{ind1,2}=folders{ind1};
	toxls{ind1,3}=cellcat(comments);
end
cd(regdir)

if exist([regdir regfile])==2
    disp('Accessing spreadsheet file - reading existing entries...')
    [exceldates,cmmnts,raw]=xlsread([regdir regfile]);
    matlabDates = datenum('30-Dec-1899') + exceldates;
    ind1=1;
    raw2=[];
    for ind2=1:size(raw,1)
        raw{ind2,1}=datestr(raw{ind2,1},0);
        flag=0;
        for ind3=1:size(toxls,1)
            if strcmp(raw(ind2,2),folders{ind3});
                flag=1;
            end
        end
        if flag==0
            raw2{ind1,1}=raw{ind2,1};
            raw2{ind1,2}=raw{ind2,2};
            raw2{ind1,3}=raw{ind2,3};
            if length(raw2{ind1,2})>maxlength
                maxlength=length(raw2{ind1,2});
            end
            ind1=ind1+1;
        end
    end
    toxls=vertcat(raw2,toxls);
end
varargout{1}=toxls;

for ind3=1:size(toxls,1)
	temp(ind3,:)=[toxls{ind3,2} blanks(maxlength-length(toxls{ind3,2}))];
end
[temp2,m,n]=unique(temp,'rows');
for ind3=1:size(toxls,1)
	temp3{n(ind3),1}=toxls{ind3,1};
	temp3{n(ind3),2}=toxls{ind3,2};
    temp3{n(ind3),3}=toxls{ind3,3};
end
disp('Accessing spreadsheet file - writing new entries...')
xlswrite(regfile,temp3);
cd(origpath)
disp('...done.')
end
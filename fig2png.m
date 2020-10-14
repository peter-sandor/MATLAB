function fig2png(varargin)
% converts all .fig files within a given directory and all its subdirectories to .png
% target dir is supplied in argument, if not, then the default is the
% current directory. There is an optional switch '0' which sets directory
% for search to be the target directory without its subdirectories.

if nargin==1
    if isstr(varargin{1})
        path=varargin{1};
        if path(end)~='\'
            path=[path '\'];
        end
        flag=1;
    else flag=varargin{1};
        path=[];
    end
elseif nargin==0
    path='.\';
    flag=1;
elseif nargin==2
    path=varargin{1};
    if path(end)~='\'
        path=[path '\'];
    end
    flag=varargin{2};
end
if exist(path)==0
    disp('Nonexistent path!');
    return
end    
if flag==1
    files=subdir(strcat(path,'*.fig'));
elseif flag==0
    files=dir(strcat(path,'*.fig'));
end
N=length(files);
disp(['Found ' num2str(N) ' ".fig" file(s)'])
for ind1=1:N
    h=openfig(files(ind1).name,'new','invisible');
    set(h,'color','white')
%     set(h,'units','normalized','outerposition',[0 0 1 1]); % make figure fullscreen
    img = getframe(h);
    imwrite(img.cdata, [files(ind1).name(1:max(strfind(files(ind1).name,'.fig'))-1) '.png']);
%     saveas(h,files(ind1).name(1:max(strfind(files(ind1).name,'.fig'))-1),'png');
    disp(strcat(num2str(ind1),'/',num2str(N)));
    close(h);
end
end
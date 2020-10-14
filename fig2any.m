function fig2any(varargin)

% fig2any               - converts .fig files in current directory to .png
% fig2any(ext)          - converts .fig files in current directory to image
%                           files specified in 'ext'
% fig2any(ext,path)     - converts .fig files in directory specified in
%                           'path' to image files specified in 'ext'

% specifiy 'epsc' as extension to save graphics in color .eps format!
if nargin==2
    path=varargin{2};
    extension=varargin{1};
elseif nargin==1
    path=[];
    extension=varargin{1};
elseif nargin==0
    path=[];
    extension='png';
end
files=dir(strcat(path,'*.fig'));
N=length(files);
disp(['Found ' num2str(N) ' ".fig" file(s)'])
for ind1=1:N
    h=openfig(strcat(path,files(ind1).name),'new','invisible');
%     print(h, '-opengl', '-depsc', '-r300', files(ind1).name(1:max(strfind(files(ind1).name,'.fig'))-1));
    saveas(h,strcat(path,files(ind1).name(1:max(strfind(files(ind1).name,'.fig'))-1)),extension);
    disp([num2str(ind1) '/' num2str(N)]);
    close(h);
end
end
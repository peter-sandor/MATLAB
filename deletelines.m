function varargout = deletelines(filenamein,N,varargin) 

fid = fopen(filenamein, 'r'); % Open source file.
header=[];
for ind1=1:N
	header=cat(1,header,str2num(fgetl(fid))); % Read/discard line.
end
buffer = fread(fid,Inf);	% Read rest of the file.
fclose(fid);
[pathstring,namein,extin] = fileparts(filenamein);
if nargin==2
    filenameout=[pathstring namein '_truncated' extin];
elseif nargin==3;
    filenameout=[pathstring varargin{1}];
end
fid = fopen(filenameout,'w');   % Open destination file.
fwrite(fid, buffer);	% Save to file.
fclose(fid);
varargout{1}=filenameout;
varargout{2}=header;
end
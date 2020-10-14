function varargout = loadbinaryimage2(filename)
[fileID,message] = fopen(filename);
if fileID==-1
    disp([filename ': ' message]);
    output=[];
else
    output = fread(fileID,[360 360],'uint8','ieee-be');
% 	output = fread(fileID,'float32' ,'ieee-be')';
    fclose(fileID);
end
varargout{1}=output;
end
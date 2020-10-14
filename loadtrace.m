function varargout = loadtrace(filename)
[fileID,message] = fopen(filename);
if fileID==-1
    disp([filename ': ' message]);
    output=[];
else
    output = fgetl(fileID);
    if output==-1;
        output=[];
    elseif ~isempty(output)
        output=str2num(output);
    end
% 	output = fread(fileID,'float32' ,'ieee-be')';
    fclose(fileID);
end
varargout{1}=output;
end
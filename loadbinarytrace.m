function output = loadbinarytrace(filename)
fileID = fopen(filename);
output = fread(fileID,[8000 1],'float32' ,'ieee-be')';
%output = fread(fileID,'float32' ,'ieee-be')';
fclose(fileID);
end
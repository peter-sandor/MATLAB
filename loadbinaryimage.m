function output = loadbinaryimage(filename)
fileID = fopen(filename);
output = fread(fileID,[1024 768],'float32' ,'ieee-be')';
%output = fread(fileID,'float32' ,'ieee-be')';
fclose(fileID);
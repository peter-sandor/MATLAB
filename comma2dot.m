function varargout = comma2dot(filename)
FID1=fopen(filename,'r');

str_to_write=[];
while ~feof(FID1)
   str1=fgets(FID1);
   index1=strfind(str1,',');
   if ~isempty(index1) && ~isnumeric(str1)
       str1(index1)='.';
   end
   str_to_write=[str_to_write str1];
end
fclose(FID1);
FID2=fopen([filename(1:end-4) '_dot.txt'],'w');
fprintf(FID2,'%s',str_to_write);
fclose(FID2);
varargout{1}=str_to_write;
end
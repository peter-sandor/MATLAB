function data = read_itx(filename)
% This function reads data from the IGOR text file (.itx).
fid = fopen(filename);
% headerrows = 9;
% format = repmat('%f ', 1, Ncols);
% data = cell2mat(textscan(fid, format, 'MultipleDelimsAsOne',true, 'Delimiter', ' ', 'headerLines', headerrows));
while ~strcmp(fgetl(fid),'BEGIN')
end
counts = map2colvec(str2num(fgetl(fid)));
N = length(counts);
fgetl(fid);
x_text = fgetl(fid);
fclose(fid);
ind_commas = strfind(x_text,',');
E_min = str2num(x_text(ind_commas(1)+1:ind_commas(2)-1));
E_max = str2num(x_text(ind_commas(2)+1:ind_commas(3)-1));
E_step = (E_max - E_min)/(N-1);
E = map2colvec(E_min:E_step:E_max);
data = [E counts];
end
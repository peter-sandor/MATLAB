function varargout = pad_name_w_zeros(filename,num_width,copy_or_rename)
% This function pads the end of the filename (supposedly an index number) with zeros to obtain a fixed width of that number ('num_width').
% Based on the string input 'copy_or_rename', it creates a copy of the file (copy_or_rename = 'copy') or replaces the existing one (copy_or_rename = 'rename').
% The function can return with the new filename
ind_space = strfind(filename,' ');
ind_uscore = strfind(filename,'_');
ind_dot = strfind(filename,'.');
ind_last = max([max(ind_space) max(ind_uscore)]);
current_width = ind_dot - ind_last - 1;
Nzeros = num_width-current_width;
str_zeros = [];
for ind1 = 1:Nzeros
    str_zeros = [str_zeros '0'];
end
newname = [filename(1:ind_last) str_zeros filename(ind_last+1:end)];
varargout{1} = newname;
if ~strcmp(newname,filename)
    if strcmp(copy_or_rename,'copy')
        copyfile(filename,newname);
    elseif strcmp(copy_or_rename,'rename')
        movefile(filename,newname);
    end
end
end
function varargout = space2uscore(filename,copy_or_rename)
% This function replaces whitespace characters in 'filename' with
% underscores. Based on the string input 'copy_or_rename', it creates a copy of the
% file (copy_or_rename = 'copy') or replaces the existing one (copy_or_rename = 'rename').
% The function can return with the new filename.
oldname = filename;
ind_space = strfind(oldname,' ');
ind_dot = strfind(oldname,'.');
if ~isempty(ind_space)
    for ind2 = 1:length(ind_space)
        newname = [oldname(1:ind_space(ind2)-1) '_' oldname(ind_space(ind2)+1:end)];
        oldname = newname;
    end
    varargout{1} = newname;
    if strcmp(copy_or_rename,'copy')
        copyfile(filename,newname);
    elseif strcmp(copy_or_rename,'rename')
        movefile(filename,newname);
    end
end
end
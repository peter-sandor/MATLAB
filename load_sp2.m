function [varargout] = load_sp2(filename)
% This function loads data from SPECS Prodigy .sp2 files. Keeping the
% header is optional.
fid = fopen(filename);
tline = fgetl(fid);
data = [];
header = [];
ind_line = 0;
while ischar(tline)
    ind_line = ind_line + 1;
    if ind_line == 1447732
        1;
    end
    tline = fgetl(fid);
    if ind_line<51
        header = [header ' ' tline];
    elseif ind_line == 51
        N = str2num(tline);
    else
        if ischar(tline) && strcmp(tline,'P2')
            break;
        end
        if ischar(tline)
            data = [data; str2num(tline)];
        end
    end
end
fclose(fid);
pic_out = reshape(data,[N(1) N(2)]).';

varargout{1} = pic_out;
varargout{2} = header;
end
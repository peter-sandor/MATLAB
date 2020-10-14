function varargout = read_itx_multiaxes(filename)
% This function reads data from the IGOR text file (.itx).
fid = fopen(filename);
% Look for place in file where data begins.
while ~strcmp(fgetl(fid),'BEGIN')
end
% Start reading data, and look for place where it ends.
ind1 = 1;
temp_text = fgetl(fid);
while ~strcmp(temp_text,'END')
    counts(ind1,:) = map2colvec(str2num(temp_text));
    temp_text = fgetl(fid);
    ind1 = ind1 + 1;
end
N_E = size(counts,2);
N_alpha = size(counts,1);

if N_alpha>1
    ind_lookfor = 2;
else
    ind_lookfor = 1;
    alpha = [];
end
% Look for lines with info on axes (energy & angle, if applicable)
for ind1 = 1:ind_lookfor
    temp_text = fgetl(fid);
    while isempty(strfind(temp_text,'X SetScale/I'))
        temp_text = fgetl(fid);
    end
    ind_quotes = strfind(temp_text,'"');
    switch temp_text(ind_quotes(1)+1:ind_quotes(2)-1)
        case 'Non-Energy Channel [mm]'
            ind_commas = strfind(temp_text,',');
            alpha_min = str2num(temp_text(ind_commas(1)+1:ind_commas(2)-1));
            alpha_max = str2num(temp_text(ind_commas(2)+1:ind_commas(3)-1));
            alpha_step = (alpha_max - alpha_min)/(N_alpha-1);
            alpha = map2colvec(alpha_min:alpha_step:alpha_max);
        case 'Non-Energy Channel [deg]'
            ind_commas = strfind(temp_text,',');
            alpha_min = str2num(temp_text(ind_commas(1)+1:ind_commas(2)-1));
            alpha_max = str2num(temp_text(ind_commas(2)+1:ind_commas(3)-1));
            alpha_step = (alpha_max - alpha_min)/(N_alpha-1);
            alpha = map2colvec(alpha_min:alpha_step:alpha_max);
        case 'Kinetic Energy [eV]'
            ind_commas = strfind(temp_text,',');
            E_min = str2num(temp_text(ind_commas(1)+1:ind_commas(2)-1));
            E_max = str2num(temp_text(ind_commas(2)+1:ind_commas(3)-1));
            E_step = (E_max - E_min)/(N_E-1);
            E = map2colvec(E_min:E_step:E_max);
        case 'cps (Intensity)'
            1; % placeholder
        otherwise
            1;
    end
end
fclose(fid);

varargout{1} = counts;
varargout{2} = E;
varargout{3} = alpha;
end
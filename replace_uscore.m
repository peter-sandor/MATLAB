function text_out = replace_uscore(text_in)
% This function takes a string input and replaces the underscores in it with
% whitespaces.

indices=strfind(text_in,'_');
text_out=text_in;
text_out(indices)=' ';
end
function [foldername,filename] = strip_path(path_in)

index1=max(vec2ind(path_in=='\'));
foldername=path_in(1:index1);
filename=path_in(index1+1:end);
end
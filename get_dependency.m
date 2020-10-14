function varargout = get_dependency(filename)
[fList,pList] = matlab.codetools.requiredFilesAndProducts(filename);
disp([pList.Name ' ' pList.Version])
for ind1 = 1:length(fList)
    disp(fList{ind1});
end
if nargout>0
    varargout{1} = fList;
    varargout{2} = pList;
end
end
function filename = getFilenameString(string)

ind_slash = regexp(string,'/');

filename = string(ind_slash(end)+1:end);
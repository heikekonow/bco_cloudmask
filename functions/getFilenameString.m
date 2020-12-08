% 	Function to extract the file name from a full path.
%
% 	Input variables:
%		- string 	path string
%
%   Output variables:
%       - filename  file name from from input string
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020

function filename = getFilenameString(string)

% Look for slashes in path string
ind_slash = regexp(string,'/');

% File name is the part after the last slash
filename = string(ind_slash(end)+1:end);
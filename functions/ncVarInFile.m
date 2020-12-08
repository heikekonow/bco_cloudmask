% 	Function to check if variable is in NetCDF file.
%
% 	Input variables:
%		- file      string with full path to NetCDF file
%       - var       string with variable name
%
%   Output variables:
%       - result    number of times the variable was found in 'file'
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020

function result = ncVarInFile(file,var)

% List all variables in file
vars = nclistvars(file);

% Count number of occurences of variables in file
result = sum(ismember(vars,var));
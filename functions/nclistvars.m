% 	Function to return variable information from NetCDF file.
%
% 	Input variables:
%		- pathtofile 	string with full path to NetCDF file
%
%   Output variables:
%       - varnames      cell array containing variable names
%       - dimids        cell array containing dimension identifiers
%       - dimnames      cell array containing dimension names
%       - dimlengths    cell array containing dimension lengths
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020


function [varnames,dimids,dimnames,dimlengths] = nclistvars(pathtofile)

% Function to list variables in netcdf file

% Open file
ncid = netcdf.open(pathtofile,'NC_NOWRITE');
% Get number of variables
[~,numvars,~,~] = netcdf.inq(ncid);

% Preallocate
varnames = cell(numvars,1);
dimids = cell(numvars,1);
dimnames = cell(1,1);
dimlengths = cell(1,1);

% Loop all entries
% netcdf.inqVar function starts counting with 0
for i = 0:(numvars - 1)
    % Get name and dimension id for variable
    [varnames{i+1,1},~,dimids{i+1,1},~] = netcdf.inqVar(ncid,i);
    % Copy dimension id
    dimIDsForVar = dimids{i+1,1};
    % Loop dimensions of this variable
    for j=1:length(dimIDsForVar)
        % Get name and length of dimensions
        [dimnames{i+1,j},dimlengths{i+1,j}] = netcdf.inqDim(ncid,dimIDsForVar(j));
    end
end

% Close file
netcdf.close(ncid)

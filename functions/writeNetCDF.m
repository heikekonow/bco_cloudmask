%% writeNetCDF
%   writeNetCDF - Function to write data to NetCDF file
%
%   Syntax:  writeNetCDF(outfile,ncVarNames,ncDims,varData,varInfo)
%
%   Inputs:
%       outfile     - String defining the full path of the file to be written
%       ncVarNames  - Cell array of variable names in workspace
%       ncDims      - Cell array of variable dimension information:
%                     dimension name(s) and dimension length(s) in one cell
%       varData     - Cell array of variable data
%       varInfo     - Cell array with variable infos: units, long names
%       globAtt     - Cell array of global attribute data
%   Set NetCDF version and compression level in beginning of script.
%
%   Outputs:
%       NetCDF file as defined in 'outfile' variable, no output variables
%
%   Example:
%
%       outfile = '~/data/superduperfile.nc';
%
%       ncVarNames = {'TIME','HEIGHT','MR'};
%
%       ncDims = {{'TIME',15889};...
%                  {'HEIGHT',412};...
%                   {'HEIGHT',412,'TIME',15889}};
%
%       varData = {1:15889;...
%                  1:412;...
%                  rand([412,15889])};
%
%       varInfo = {...
%                   'time','seconds since 1970-01-01 00:00:00 UTC','TIME';
%                   'height','m','HEIGHT';
%                   'mixratio','g/kg','MR'};
%
%       globAtt = {{'contact','heike.konow@uni-hamburg.de'};...
%                   {'mission','NAWDEX'};...
%                   {'flight_number',10}};
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: 
%
%   Author: Dr. Heike Konow
%   Meteorological Institute, Hamburg University
%   email address: heike.konow@uni-hamburg.de
%   Website: http://www.mi.uni-hamburg.de/
%   August 2016; Last revision: January 2017

%%
function writeNetCDF(outfile,ncVarNames,ncDims,varData,varInfo,globAtt)

%------------- BEGIN CODE --------------

%%%%% Input %%%%%%%%%%%
% Define netcdf format
ncformat = 'netcdf4';
% Define netcdf compression (value between 0 and 9, with 9 beeing strongest compression) 
deflatelevel = 5;
%%%%%%%%%%%%%%%%%%%%%%%

% Delete netcdf file if it already exists
if exist(outfile,'file')
    delete(outfile)
end
%%

% Loop all variables
for i=1:length(ncVarNames)
    
    % Display progress
    disp(['Writing (' num2str(i) '/' num2str(length(ncVarNames)) '):' ncVarNames{i}])
    
    % Find index of variable in third or fourth column of info matrix
    if size(varInfo,2)==4
        ind = strcmp(varInfo(:,4),ncVarNames{i});
    else
        ind = strcmp(varInfo(:,3),ncVarNames{i});
    end
    
    % check if variable is time and in fact sdn (i.e. < 800000), in this
    % case: convert to unix time
    if strcmp(ncVarNames{i},'time') && varData{i}(1)<800000 
        data = varData{i};
        data = sdn2unixtime(data);
        
    % Else, if variable is time, and alread unix time, round to integer and
    % check afterwards
    elseif strcmp(ncVarNames{i},'time') && varData{i}(1)>800000 
        data = checkAndRoundTime(varData{i});
    else
        data = varData{i};
%         data = single(varData{i});
    end
    
    % If variable is scalar, -> only one dimension
    if length(ncDims{i})==1
        % Create netcdf variable
        nccreate(outfile,ncVarNames{i},'Format',ncformat,'DeflateLevel',deflatelevel)
        
    % If variable is vector or matrix
    else
        % Create netcdf variable
        
%         if ~strcmp(ncVarNames{i},'time') && varData{i}(1)<800000 
%             nccreate(outfile,ncVarNames{i},'Dimensions',[ncDims{i}],'Datatype','single',...
%                 'Format',ncformat,'DeflateLevel',deflatelevel)
%         else % if variable is 'time'
            nccreate(outfile,ncVarNames{i},'Dimensions',[ncDims{i}],'Datatype','double',...
                'Format',ncformat,'DeflateLevel',deflatelevel)
%         end
    end
    
    % Write data to file
    ncwrite(outfile,ncVarNames{i},data)
    
    % Write unit attribute
    ncwriteatt(outfile,ncVarNames{i},'units',varInfo{ind,2})
    
    % Write long_name attribute
    ncwriteatt(outfile,ncVarNames{i},'long_name',varInfo{ind,1})
    
    % Clean up
    clear data
end

% Write global attributes to netcdf file
for i=1:length(globAtt)
    ncwriteatt(outfile,'/',globAtt{i}{1},globAtt{i}{2})
end

%------------- END OF CODE --------------
end

function timeRounded = checkAndRoundTime(time)
    timeRounded = round(time);
    
    if numel(time)~=numel(unique(timeRounded))
        error('rounded time is funny')
    end
end

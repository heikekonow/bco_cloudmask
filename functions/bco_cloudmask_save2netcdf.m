% 	Code that looks for and concatenates bco data (wind and radar).
% 	Input variables:
%		- start_date	string with first date to process (yyyymmdd format)
%		- end_date		string with last date to process (yyyymmdd format)
%		- radarname		string with radar name (MBR or KATRIN)
%		- radarrange	string with radar height range (e.g. '155m-18m')
%       - vers          string with version number (e.g. v0.3)
%       - newextra      set to true if additional data should be processed and saved,
%                       additional variables are reflectivity, cloud mask,
%                       closed cloud mask, wind speed, LDR, Doppler velocity,
%                       spectral width
%		- instrument 	instrument name for output netcdf globa attribute
%		- outfolder     output directory
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020



function bco_cloudmask_save2netcdf(start_date, end_date, radarname, radarrange, vers, newextra, instrument, outfolder, tmppath)

% Set paths to temporary file and output files
filepath = [tmppath 'Z_' radarname '_' radarrange '_' start_date '-' end_date '_closed_concomp.mat'];
outfile = [outfolder 'cloudObjectMask_' radarname '_' radarrange '_' start_date '-' end_date '_' vers '.nc'];
outfile_2 = [outfolder 'cloudObjectMask_' radarname '_' radarrange '_' start_date '-' end_date '_extradata_' vers '.nc'];

% Display
disp('>>>>>>>>>> <<<<<<<<<<')
disp('Output file:')
disp(outfile)
disp('>>>>>>>>>> <<<<<<<<<<')

% List of variable names in temporary file
varnames2read = {...%'Z',
                'cloudBase','cloudDepth','cloudEndTime','cloudLength',...%'cloudMask','cloudMask_closed',
                'cloudStartTime','cloudTop',...%'windSpeed',...
                'height','numMask','numberOfClouds','time','status','wind_missing'};

% Display
disp('Reading saved data from mat file')

% Load data
load(filepath,varnames2read{:})

% List of variable names in netcdf file
ncVarNames = {...%'Zf',
              'cloudBase','cloudDepth','cloudEndTime','cloudLength',...%'cloudMask','cloudMask_closed',
              'cloudStartTime','cloudTop',...%'windSpeed',...
              'height','numMask','N_clouds','time','status', 'wind_missing'};

% Concatenate variable data
varData = {...%Z,
           cloudBase{1}, cloudDepth{1}, cloudEndTime{1}, cloudLength{1}, ...%cloudMask,cloudMask_closed,
           cloudStartTime{1}, cloudTop{1},... % windSpeed{1},...
           height{1}, numMask, 1:numberOfClouds ,time{1}, status, wind_missing};

% List variable dimensions
ncDims = {...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'height',length(height{1})};...
          {'height',length(height{1}),'time',length(time{1})};...
          {'N_clouds',numberOfClouds};...
          {'time',length(time{1})};...
          {'time',length(time{1})};...
          {'time', length(time{1})},...
          };

% List variable units
units = {'m';
         'm';
         'seconds since 1970-1-1 0:00:00 UTC';
         'm';
         'seconds since 1970-1-1 0:00:00 UTC';
         'm';
         'm';
         ' ';
         ' ';
         'seconds since 1970-1-1 0:00:00 UTC';
         ' ';
         ' '};

% List variable long names
long_name = {...
             'Cloud base height';
             'Cloud depth';
             'End of cloud';
             'Cloud length';
             'Beginning of cloud';
             'Cloud top height';
             'Height';
             'Cloud object ID';
             'Total number of clouds';
             'Time';
             'Radar On (1), radar Off (0), no data file (2), radar scanning (3)';
             'Wind measured (0), wind missing (1)'};

% Concatenate all information about variable that will be written to netcdf
varInfo = [long_name, units, varnames2read', ncVarNames'];

% List of global attributes
globAtt = {{'information','Cloud mask applied to BCO radar'};...
        {'contact','heike.konow@mpimet.mpg.de'};...
        {'location','The Barbados Cloud Observatory, Deebles Point, Barbados, West Indies'};...
        {'institution','Max Planck Institute for Meteorology, Hamburg'};...
        {'instrument','MBR2 cloud radar'}};

% Write data to netcdf file
writeNetCDF(outfile,ncVarNames,ncDims,varData,varInfo,globAtt)

% If an additional file with extra data should be written
if newextra

    % List of variable names in temporary file
    varnames2read_2 = {'Z','cloudMask','cloudMask_closed','windSpeed',...
    'LDR', 'VEL', 'RMS'};

    % Load data
    load(filepath,varnames2read_2{:})

    % List of variable names in netcdf file
    ncVarNames_2 = {'Zf','cloudMask','cloudMask_closed','windSpeed',...
    'LDR', 'VEL', 'RMS'};

    % Concatenate variable data
    varData_2 = {Z,cloudMask,cloudMask_closed, windSpeed{1}, LDR, VEL, RMS};

    % List variable dimensions
    ncDims_2 = {{'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})};...
           {'height',length(height{1}),'time',length(time{1})}
            };

    % List variable units
    units_2 = {'dBZ'; ''; ''; 'm/s'; 'dbZ'; 'm/s'; 'm/s'};

    % List variable long names
    long_name_2 = {'Filtered and Mie corrected Radar Reflectivity of all Hydrometeors';
                'Cloud mask';
                'Cloud mask after morphological closing';
                'Wind speed';
                'Linear De-Polarization Ratio LDR of all Hydrometeors';
                'Doppler Velocity VEL of all Hydrometeors';
                'Peak Width RMS of all Hydrometeors'};
    % Concatenate all information about variable that will be written to netcdf
    varInfo_2 = [long_name_2, units_2, varnames2read_2', ncVarNames_2'];

    % List of global attributes
    globAtt_2 = {{'information','Additional data for BCO radar cloud mask'};...
            {'contact','heike.konow@mpimet.mpg.de'};...
            {'location','The Barbados Cloud Observatory, Deebles Point, Barbados, West Indies'};...
            {'institution','Max Planck Institute for Meteorology, Hamburg'};...
            {'instrument',[instrument ' cloud radar']}};

    % Write data to netcdf file
    writeNetCDF(outfile_2,ncVarNames_2,ncDims_2,varData_2,varInfo_2,globAtt_2)
end

disp('Data saved to netcdf')

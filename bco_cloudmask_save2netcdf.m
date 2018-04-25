function bco_cloudmask_save2netcdf(start_date, end_date, radarname, radarrange, version)

% clear; close all

filepath = ['/scratch/local1/m300512/bco_concat/Z_' radarname '_' radarrange '_' start_date '-' end_date '_closed_concomp.mat'];
outfile = ['~/bco_cloudmask/cloudObjectMask_' radarname '_' radarrange '_' start_date '-' end_date '_' version '.nc'];
outfile_2 = ['~/bco_cloudmask/cloudObjectMask_' radarname '_' radarrange '_' start_date '-' end_date '_extradata_v0.1.nc'];

varnames2read = {...%'Z',
                'cloudBase','cloudDepth','cloudEndTime','cloudLength',...%'cloudMask','cloudMask_closed',
                'cloudStartTime','cloudTop',...%'windSpeed',...
                'height','numMask','numberOfClouds','time','status','wind_missing'};

varnames2read_2 = {'Z','cloudMask','cloudMask_closed','windSpeed'};

disp('Reading saved data from mat file')
load(filepath,varnames2read{:})
load(filepath,varnames2read_2{:})

ncVarNames = {...%'Zf',
              'cloudBase','cloudDepth','cloudEndTime','cloudLength',...%'cloudMask','cloudMask_closed',
              'cloudStartTime','cloudTop',...%'windSpeed',...
              'height','numMask','N_clouds','time','status', 'wind_missing'};

ncVarNames_2 = {'Zf','cloudMask','cloudMask_closed','windSpeed'};

varData = {...%Z,
           cloudBase{1}, cloudDepth{1}, cloudEndTime{1}, cloudLength{1}, ...%cloudMask,cloudMask_closed,
           cloudStartTime{1}, cloudTop{1},... % windSpeed{1},...
           height{1}, numMask, 1:numberOfClouds ,time{1}, status, wind_missing};

varData_2 = {Z,cloudMask,cloudMask_closed, windSpeed{1}};

ncDims = {...%'height',length(height{1}),'time',length(time{1})};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
%          {'height',length(height{1}),'time',length(time{1})};...
%          {'height',length(height{1}),'time',length(time{1})};...
          {'N_clouds',numberOfClouds};...
          {'N_clouds',numberOfClouds};...
%          {'height',length(height{1}),'time',length(time{1})};...
          {'height',length(height{1})};...
          {'height',length(height{1}),'time',length(time{1})};...
          {'N_clouds',numberOfClouds};...
          {'time',length(time{1})};...
          {'time',length(time{1})};...
          {'time', length(time{1})},...
          };


ncDims_2 = {{'height',length(height{1}),'time',length(time{1})};...
       {'height',length(height{1}),'time',length(time{1})};...
       {'height',length(height{1}),'time',length(time{1})};...
       {'height',length(height{1}),'time',length(time{1})};...
        };

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

units_2 = {'dBZ'; ''; ''; 'm/s'};

long_name = {...%'Filtered and Mie corrected Radar Reflectivity of all Hydrometeors';
             'Cloud base height';
             'Cloud depth';
             'End of cloud';
             'Cloud length';
%             'Cloud mask';
%             'Cloud mask after morphological closing';
             'Beginning of cloud';
             'Cloud top height';
%             'Wind speed';
             'Height';
             'Cloud object ID';
             'Total number of clouds';
             'Time';
             'Radar On (1), radar Off (0), no data file (2), radar scanning (3)';
             'Wind measured (0), wind missing (1)'};

long_name_2 = {'Filtered and Mie corrected Radar Reflectivity of all Hydrometeors';
            'Cloud mask';
            'Cloud mask after morphological closing';
            'Wind speed';};

varInfo = [long_name, units, varnames2read', ncVarNames'];
varInfo_2 = [long_name_2, units_2, varnames2read_2', ncVarNames_2'];

globAtt = {{'information','Cloud mask applied to BCO radar'};...
        {'contact','heike.konow@mpimet.mpg.de'};...
        {'location','The Barbados Cloud Observatory, Deebles Point, Barbados, West Indies'};...
        {'institution','Max Planck Institute for Meteorology, Hamburg'};...
        {'instrument','MBR2 cloud radar'}};


globAtt_2 = {{'information','Additional data for BCO radar cloud mask'};...
        {'contact','heike.konow@mpimet.mpg.de'};...
        {'location','The Barbados Cloud Observatory, Deebles Point, Barbados, West Indies'};...
        {'institution','Max Planck Institute for Meteorology, Hamburg'};...
        {'instrument','MBR2 cloud radar'}};

writeNetCDF(outfile,ncVarNames,ncDims,varData,varInfo,globAtt)
writeNetCDF(outfile_2,ncVarNames_2,ncDims_2,varData_2,varInfo_2,globAtt_2)

disp('Data saved to netcdf')

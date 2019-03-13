
% Time
tic

% Clean up
clear; close all

% disp('hello')
% pause;

% Set path to radar data
path = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/W-Band_Radar/';
% Set path to output files
outpath = '/scratch/local1/m300512/bco_cloudmask_wband/';
% outpath = '/pool/OBS/ACPC/MBR2/cloudmask/bco_object_cloudmask/cloudObjectMask';

% List monthly folders
monthfolders = listFiles([path '20*'], 'full');

% Set radar names to work on
radarname = {'W-Band'};
% Set version for output nc file
version = 'v0.3';
% Write new version of additional data file?
newextra = true;
% Set threshold for reflectivity
dbz_threshold = -500;
% Set value for radar range
radarrange = 'na'


%% Prepare dates %%%%%%%%%%%%%%%%%%%%%%%%

% Loop months
for i=1:length(monthfolders) %%% !!!!!!!!!!!!!!!

    % List files in monthly folder
    radarfiles = listFiles([monthfolders{i} '/*.nc']);

    % Analyse file name to find double underscores to extract parts of file names
    b{i} = cellfun(@(x) regexp(x, '__'), radarfiles, 'uni', false);

    % Extract dates from file names
    dates{i} = cell2mat(cellfun(@(x,y) x(y(end)+2:end-3), radarfiles, b{i}, 'uni', false));

    % Display info
    disp(['>>>>>>>> Processing ' radarname{1} ' radar for ' monthfolders{i} ' <<<<<<<<'])

    % Get first and last date of data
    %%% !!!!!!!!!!!!!!!!!!!!!!!
    start_date = num2str(min(str2num(dates{i})));
    % start_date = dates{i}(19, :);
    %%% !!!!!!!!!!!!!!!!!!!!!!!
    end_date = num2str(max(str2num(dates{i})));

    %% Actual processing %%%%%%%%%%%%%%%%%%%%%%%%

    % Generate output file names to test if files exist already
    outfile = [outpath '_' radarname{1} '_' start_date '-' end_date '_' version '.nc'];
    outfile_2 = [outpath '_' radarname{1} '_' start_date '-' end_date '_extradata_' version '.nc'];

    % Copy variable
    datafiles = radarfiles;

    % Only do processing if outfile doesn't already exist and
    % ignore if string 'deg' is part of file name
    if ~exist(outfile, 'file')

        % Concatenate data
        bco_cloudmask_concatData([monthfolders{i} '/'], datafiles, radarname{1}, radarrange, start_date, end_date, dbz_threshold)

        % Generate cloud mask
        bco_cloudmask_mask(radarname{1}, radarrange, start_date, end_date)

        % Caclulate cloud parameter
        bco_cloudmask_param(start_date, end_date, radarname{1}, radarrange)

        % Save data to netcdf
        bco_cloudmask_save2netcdf(start_date, end_date, radarname{1}, radarrange, version, newextra, radarname{1}, outpath)

    % If new extra data should be processed and the corresponding file doesn't exist already
    elseif newextra && ~exist(outfile_2, 'file') && ~isempty(datafiles) %~contains(start_date, 'deg')
        bco_cloudmask_save2netcdf(start_date, end_date, radarname{1}, radarrange, version, newextra, radarname{1}, outpath)
    end


end

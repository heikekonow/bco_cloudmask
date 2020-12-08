% 	Code that starts cloud mask processing.
% 	Set data paths and processing features in the beginning of this file
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020

% Time
tic

% Clean up
clear; close all

%% Set parameters here %%%%%%%%%%%%%%%%%%%%%%%%
% % Set path to radar data
% path = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/Ka-Band/MBR2/10s/';
% % Set path to output files
% outpath = '/pool/OBS/ACPC/MBR2/cloudmask/bco_object_cloudmask/';
% Set radar names to work on
radarname = {'MBR'};
% radarname = {'MBR', 'KATRIN'};
% Set version for output nc file
version = 'v0.3';
% Write new version of additional data file?
newextra = true;
% Set threshold for reflectivity
dbz_threshold = -50;
% Set zenith angle, everything deviating from this will be removed as scanning
zenithAngle = 0;
% Set minimum cloud object size
minSize = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

configstruct = config;
path = configstruct.path;
outpath = configstruct.outpath;
outpathtmp = configstruct.outpathtmp;
tmppath = configstruct.tmppath;
windpath = configstruct.windpath;

%% Prepare dates %%%%%%%%%%%%%%%%%%%%%%%%

% Loop radars
for i=1:length(radarname)

    % List available month folders
    monthfolders = listFiles(path, 'full');

    % Look for nc files in folder list and remove
    ind_nc = contains(monthfolders, '.nc');
    monthfolders(ind_nc) = [];

    % Loop month folders
    for m=1:length(monthfolders)
        % List netcdf files in month folder
        files{m} = listFiles([monthfolders{m} '/'], 'full');
    end

    % Concatenate all files into one cell array
    files = vertcat(files{:});

    % Identify empty folders and delete those entries
    ind = cell2mat(cellfun(@isempty, files, 'uni', 0));
    files(ind) = [];

    % Analyse file name to find double underscores to extract parts of file names
    a{i} = cellfun(@(x) regexp(x, '__'), files, 'uni', false);
    % Extract height range from file names
    heightrange{i} = cellfun(@(x,y) x(y(4)+2:y(5)-1), files, a{i}, 'uni', false);
    % List unique height ranges
    unique_height{i} = unique(heightrange{i});

    dates{i} = cell2mat(cellfun(@(x,y) x(y(end)+2:end-3), files, a{i}, 'uni', false));

    % Loop unique height ranges
    for j=1:length(unique_height{i})

        % Rename variable
        radarrange = unique_height{i}{j};

        ind_heights = contains(files, radarrange);

        % List all files that match radar name and height range
        radarfiles{i,j} = files(ind_heights);
        % radarfiles{i,j} = listFiles([path 'MMCR*' radarname{i} '*' radarrange '*.nc']);

        % Analyse file name to find double underscores to extract parts of file names
        b{i,j} = cellfun(@(x) regexp(x, '__'), radarfiles{i,j}, 'uni', false);

        % Extract dates from file names
        % dates{i,j} = cell2mat(cellfun(@(x,y) x(y(5)+2:end-3), radarfiles{i,j}, b{i,j}, 'uni', false));
        % 20190206: change index of double underscores to last one (date should always be last
        %           in filenames)
        dates{i,j} = cell2mat(cellfun(@(x,y) x(y(end)+2:end-3), radarfiles{i,j}, b{i,j}, 'uni', false));

        % Display info
        disp(['>>>>>>>> Processing ' radarname{i} ' radar for ' radarrange ' range <<<<<<<<'])

        % Generate year vector strings by adding '20' in the beginning
        years = [repmat('20', size(dates{i,j},1), 1)  dates{i,j}(:,1:2)];

        % List unique years
        u_years = unique(cellstr(years));

        % Loop years
        for k=1:length(u_years)

            % Index for years matching year of current loop
            ind_years = strncmp(cellstr((u_years(k))), cellstr(years), 4);

            % List of files to read
            datafiles = radarfiles{i,j}(ind_years);

            % Generate flag to indicate not zenith pointing data
            notzenith = cellfun(@(x) contains(x, 'deg'), datafiles);

            % Delete datafiles with not zenith pointing data
            datafiles(notzenith) = [];

            % List of dates of files to read
            dates_use = dates{i,j}(ind_years,:);

            % 20190206: changed calcuation of first and last date to min and max
            %           to not be dependend on order of files on file system
            start_date = num2str(min(str2num(dates_use)));
            end_date = num2str(max(str2num(dates_use)));
            % OLD: First and last date of files to read (needed for file naming)
            % start_date = dates_use(1,:);
            % end_date = dates_use(end,:);
            % end_date = dates_use(5,:); % use for quick debugging

            %% Actual processing %%%%%%%%%%%%%%%%%%%%%%%%

            % Generate output file names to test if files exist already
            outfile = [outpath 'cloudObjectMask_' radarname{i} '_' radarrange '_' start_date '-' end_date '_' version '.nc'];
            outfile_2 = [outpath 'cloudObjectMask_' radarname{i} '_' radarrange '_' start_date '-' end_date '_extradata_' version '.nc'];

            % Only do processing if outfile doesn't already exist and
            % ignore if string 'deg' is part of file name
            if ~exist(outfile, 'file') && ~isempty(datafiles) %~contains(start_date, 'deg')
                % Concatenate data
                % bco_cloudmask_concatData(path, datafiles, radarname{i}, radarrange, start_date, end_date, dbz_threshold)
                bco_cloudmask_concatData(datafiles, radarname{i}, radarrange, start_date, end_date, dbz_threshold, zenithAngle, ...
                                            outpathtmp, tmppath, windpath)

                % Generate cloud mask
                bco_cloudmask_mask(radarname{i}, radarrange, start_date, end_date, minSize, tmppath)

                % Caclulate cloud parameter
                bco_cloudmask_param(start_date, end_date, radarname{i}, radarrange, tmppath)

                % Save data to netcdf
                bco_cloudmask_save2netcdf(start_date, end_date, radarname{i}, radarrange, version, newextra, radarname{i}, outpath, tmppath)

            % If new extra data should be processed and the corresponding file doesn't exist already
            elseif newextra && ~exist(outfile_2, 'file') && ~isempty(datafiles) %~contains(start_date, 'deg')
                bco_cloudmask_save2netcdf(start_date, end_date, radarname{i}, radarrange, version, newextra, radarname{i}, outpath)
            end

        end
    end
end

toc

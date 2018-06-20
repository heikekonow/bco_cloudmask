% Clean up
clear; close all

% Set path to radar data
path = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/Version_2/';
% Set radar names to work on
radarname = {'MBR', 'KATRIN'};
% Set version for output nc file
version = 'v0.3';
% Write new version of additional data file?
newextra = true;

%% Prepare dates %%%%%%%%%%%%%%%%%%%%%%%%

% Loop radars
for i=1:length(radarname)

    % List all files matching radar name
    files = listFiles([path '*' radarname{i} '*.nc']);
    % Analyse file name to find double underscores to extract parts of file names
    a{i} = cellfun(@(x) regexp(x, '__'), files, 'uni', false);
    % Extract height range from file names
    heightrange{i} = cellfun(@(x,y) x(y(4)+2:y(5)-1), files, a{i}, 'uni', false);
    % List unique height ranges
    unique_height{i} = unique(heightrange{i});

    % Loop unique height ranges
    for j=1:length(unique_height{i})

        % List all files that match radar name and height range
        radarfiles{i,j} = listFiles([path '*' radarname{i} '*' unique_height{i}{j} '*.nc']);

        % Analyse file name to find double underscores to extract parts of file names
        b{i,j} = cellfun(@(x) regexp(x, '__'), radarfiles{i,j}, 'uni', false);

        % Extract dates from file names
        dates{i,j} = cell2mat(cellfun(@(x,y) x(y(5)+2:end-3), radarfiles{i,j}, b{i,j}, 'uni', false));

        % Display info
        disp(['>>>>>>>> Processing ' radarname{i} ' radar for ' unique_height{i}{j} ' range <<<<<<<<'])

        % Generate year vector strings by adding '20' in the beginning
        years = [repmat('20', size(dates{i,j},1), 1)  dates{i,j}(:,1:2)];

        % List unique years
        u_years = unique(cellstr(years));

        % Loop years
        for k=1:length(u_years)

            % Index for years matchin year of current loop
            ind_years = strncmp(cellstr((u_years(k))), cellstr(years), 4);

            % List of files to read
            datafiles = radarfiles{i,j}(ind_years);

            % List of dates of files to read
            dates_use = dates{i,j}(ind_years,:);

            % First and last date of files to read (neede for file naming)
            start_date = dates_use(1,:);
            % end_date = dates_use(5,:); % use for quick debugging
            end_date = dates_use(end,:);

            %% Actual processing %%%%%%%%%%%%%%%%%%%%%%%%

            % Concatenate data
            bco_cloudmask_concatData(path, datafiles, radarname{i}, unique_height{i}{j}, start_date, end_date)

            % Generate cloud mask
            bco_cloudmask_mask(radarname{i}, unique_height{i}{j}, start_date, end_date)

            % Caclulate cloud parameter
            bco_cloudmask_param(start_date, end_date, radarname{i}, unique_height{i}{j})

            % Save data to netcdf
            bco_cloudmask_save2netcdf(start_date, end_date, radarname{i}, unique_height{i}{j}, version, newextra, radarname{i})
        end
    end
end

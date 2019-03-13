function bco_cloudmask_concatData(filepath, radarfiles, radarname, radarrange, start_date, end_date, dbz_threshold)

% Variables for debugging
a = 0; b = 0; c = 0; d = 0; f = 0; g = 0; h = 0; j = 0; k = 0; l = 0; m = 0; n = 0; o = 0; p = 0; q = 0; r = 0; s = 0; t = 0; u = 0; v = 0; w = 0; x = 0; y = 0; z = 0;

% Set upper limit for analysis
% >>> remove later?
height_limit = 13500;

% Get first and last day of dataset
% start_date = dates(1,:);
% end_date = dates(end,:);

% Set paths to data
outpath = '/scratch/local1/m300512/bco_concat/';
tmppath = '/scratch/local1/m300512/data_tmp/';
windpath = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/I_Meteorology_2m/';

% If tmppath doesn't exist, create
if ~exist(tmppath,'dir')
	mkdir(tmppath)
end

% Combine 2-digit year with '20' to get four-digit year
if str2num(start_date)<20000000
	year = [num2str(20) start_date(1:2)];
else
	year = start_date(1:4);
end

% Inform user
disp(['Starting processing for ' year])
disp('Concatenating')

% Concatenate folder path and filenames to get list of files
files = cellstr([repmat(char(filepath), length(radarfiles), 1)   char(radarfiles)]);
% Get list of folders with wind data (there's a folder for each month)
wind_folders = listFiles([windpath year '*'],'full');

% Loop all monthly wind folders and list nc files in them to get a list of
% all wind data files and their paths
for f=1:length(wind_folders)
	wind_filepaths{f} = listFiles([wind_folders{f} '/' '*.nc*'],'full');
end

% Convert file paths from cell to list
wind_files = vertcat(wind_filepaths{:});

% Generate vector of days between start_date and end_date
dayvector = datenum(start_date, 'yyyymmdd'):datenum(end_date, 'yyyymmdd');

% Preallocate arrays
Zcell = cell(length(dayvector),1);
VELcell = cell(length(dayvector),1);
RMScell = cell(length(dayvector),1);
LDRcell = cell(length(dayvector),1);
t = cell(length(dayvector),1);
h = cell(length(dayvector),1);
date = cell(length(dayvector),1);
status = cell(length(dayvector),1);
wind = cell(length(dayvector),1);
elv = cell(length(dayvector),1);
redodims = zeros(1, length(dayvector));

%% Read data

% Loop all days
for i=1:length(dayvector)
    % Look for files from current date
	foundfiles = strfind(files, datestr(dayvector(i), 'yymmdd'));
    % Get index of found files
	ind_foundfiles = cellfun(@(x) ~isempty(x), foundfiles);

    % Look for wind data files from current date
	wind_foundfiles = strfind(wind_files, datestr(dayvector(i), 'yymmdd'));
    % Get index of found wind files
	ind_wind_foundfiles = cellfun(@(x) ~isempty(x), wind_foundfiles);

	% Make time vector with 10 second intervals (needed first time in if clause)
	tgoal = dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59);

    % If radar files exist
	if [sum(ind_foundfiles)==1]

        %
        % %%% RADAR DATA FILE HANDLING %%%
        %

        % Copy name of file to variable
		filename = files{ind_foundfiles};
        % Print out file name
		disp(filename)

        % %%%
        % Check if the original file on share is zipped or unzipped. If
        % it's zipped, unzip to my own folder on scratch, if directly read
        % the .nc from share.
        % %%%
        % List files with file name in temporary folder (mine on scratch),
        % ignoring the last four characters (.bz2)
		lookforunzipped = listFiles([tmppath filename(36:end-4)]);
		% Check if files are zipped and unzip
        % If no file without ending .bz2 exists AND there is no file with
        % the ending .nc -> unzip file
        if isempty(lookforunzipped) && ~strcmp(filename(end-2:end), '.nc')
            % Display
			disp('unzipping')

            % Unzip file
			eval(['! bunzip2 -c ' filename ' > ' tmppath filename(36:end-4)])

			% Rename file in list (remove extension .bz2) add path to
			% folder on scratch
			files{ind_foundfiles} = [tmppath filename(36:end-4)];

            % Check if file can be found
			if isempty(listFiles(files{ind_foundfiles}))
				error(['Error: Something''s gone wrong with renaming of unzipped file: ' files{i}])
            end

        % If the file that was found on share already has the ending .nc
		elseif strcmp(filename(end-2:end), '.nc')
            % do nothing...
			files{ind_foundfiles} = files{ind_foundfiles};

        % if the data on the share is zipped but an unzipped version exists
        % on scratch
        else
			files{ind_foundfiles} = [tmppath filename(36:end-4)];
        end

        %
        % %%% WIND DATA FILE HANDLING %%%
        %

        % Same for wind data: look for files and unzip if needed

        % If a file with wind data exists
		if [sum(ind_wind_foundfiles)==1]

            % Copy file name to variable
			wind_filename = wind_files{ind_wind_foundfiles};

            % Check if current wind data file has ending .nc
	        lookforunzipped = strcmp(wind_filename(end-1:end), 'nc');

            % If there is no file with ending .nc AND there is no wind file
            % in temporary folder (mine on scratch)
	        if ~lookforunzipped && ...
                    isempty(listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']))

                % Display
                disp('unzipping wind')

                % ~isempty(strfind(wind_filename, 'bz')) [old version new version below recommended by matlab]
				% If the file name has ending .bz or .bz2
                if contains(wind_filename, 'bz')

                    % Unzip
					eval(['! bunzip2 -c ' wind_filename ' > ' tmppath wind_filename(70:end-4)])

					% Rename file in list (remove extension .bz2)
					wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end-4)];
                else % If the file ending is supposedly .nc
					eval([' ! cp ' wind_filename ' ' tmppath])

					% Rename file in list (set new path)
					wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end)];
				end

				% Check if file can be found
                if isempty(listFiles(wind_files{ind_wind_foundfiles}))
					error(['Error: Something''s gone wrong with renaming of unzipped file: ' files{i}])
                end

            % If there is a file with ending .nc (implicit) AND there is no
            % wind file in temporary folder (mine on scratch)
            elseif isempty(listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']))

                % Display
                disp('copying wind')

                % Copy file to temporary folder
				eval([' ! cp ' wind_filename ' ' tmppath])

				% Rename file in list (set new path)
				wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end)];

            % If there is a wind file in temporary folder
            else
				% Look for file name in temporary folder
				wind_out = listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']);
                % Rename file in list (set new path)
				wind_files{ind_wind_foundfiles} = [tmppath wind_out{1}];
            end

            %
            % %%% READ WIND DATA %%%
            %

			% Read wind velocity
	        windread = ncread(wind_files{ind_wind_foundfiles}, 'VEL');
            % Read time of wind data
	        twindread = unixtime2sdn(ncread(wind_files{ind_wind_foundfiles}, 'time'));
		else
			windread = nan(length(tgoal),1);
			twindread = tgoal;
        end

        %
        % %%% READ RADAR DATA %%%
        %

		% If variable Zf exists in file
		if ncVarInFile(files{ind_foundfiles},'Zf')
			Zread= ncread(files{ind_foundfiles},'Zf');

        % Else, check if variable Z exists in file
		elseif ncVarInFile(files{ind_foundfiles},'Z')
			Zread= ncread(files{ind_foundfiles},'Z');

		elseif strcmp(radarname, 'W-Band') && ncVarInFile(files{ind_foundfiles},'Ze')
			Zread= ncread(files{ind_foundfiles},'Ze');
        % Else, return an error since the reflectivity data is missing
        else
			error(['Variable Z or Zf not found in file' files{ind_foundfiles}])
        end

        % Look if additional variables exist, otherwise fill with nans
        if ncVarInFile(files{ind_foundfiles},'VEL')
            VELread = ncread(files{ind_foundfiles}, 'VEL');
		elseif strcmp(radarname, 'W-Band') && ncVarInFile(files{ind_foundfiles},'MeanVel')
            VELread = ncread(files{ind_foundfiles}, 'MeanVel');
        else
            VELread = nan(size(Zread));
        end
        if ncVarInFile(files{ind_foundfiles},'RMS')
            RMSread = ncread(files{ind_foundfiles}, 'RMS');
		elseif strcmp(radarname, 'W-Band') && ncVarInFile(files{ind_foundfiles},'SpecWidth')
            RMSread = ncread(files{ind_foundfiles}, 'SpecWidth');
        else
            RMSread = nan(size(Zread));
        end
        if ncVarInFile(files{ind_foundfiles},'LDR')
            LDRread = ncread(files{ind_foundfiles}, 'LDR');
        else
            LDRread = nan(size(Zread));
        end

        % Read time data
		tread = unixtime2sdn(ncread(files{ind_foundfiles}, 'time'));
        % Read radar status
		if ncVarInFile(files{ind_foundfiles},'status')
			statusread = ncread(files{ind_foundfiles}, 'status');
		elseif strcmp(radarname, 'W-Band') && ncVarInFile(files{ind_foundfiles},'Status')
			statusread = ncread(files{ind_foundfiles}, 'Status');
		end
        % Read radar elevation
		if ncVarInFile(files{ind_foundfiles},'elv')
        	elvread = ncread(files{ind_foundfiles}, 'elv');
		elseif strcmp(radarname, 'W-Band') && ncVarInFile(files{ind_foundfiles},'Inc_El')
			elvread = ones(length(statusread), 1) .* 90;
			% elvread = ncread(files{ind_foundfiles}, 'Inc_El');
		end
        % Read range data
		h{i} = ncread(files{ind_foundfiles},'range');

        %
        % %%% UNIFY DATA %%%
        %

        % Get reflectivity data onto uniform grid (function defined at the bottom)
        Zcell{i} = fillData(tgoal, tread, Zread);
        % Do the same for velocity, rms, ldr
        VELcell{i} = fillData(tgoal, tread, VELread);
        RMScell{i} = fillData(tgoal, tread, RMSread);
        LDRcell{i} = fillData(tgoal, tread, LDRread);
        % Get status data onto uniform grid (function defined at the bottom)
		status{i} = fillData(tgoal, tread, statusread);
        % Get wind data onto uniform grid (function defined at the bottom)
        wind{i} = fillData(tgoal, twindread, windread);
        % Get elevation data onto uniform grid (function defined at the bottom)
		elv{i} = fillData(tgoal, tread, elvread);

		% If variable elv{i} is not only filled with nans
		if sum(isnan(elv{i}))~=length(elv{i})
			% Interpolate elevation data to fill nans resulting from regridding
			elv{i} = naninterp(elv{i});

		% If variable elv{i} only consists of nans
		else
			% Set elevation to 90 deg. It's only an assumption but thats the best we know
			elv{i} = repmat(90, 1, length(elv{i}));
		end
		% Replace elevation data again with nan when the radar was not operating (status = 0)
		elv{i}(status{i}==0) = nan;

        % Copy output time array to overall time cell array
        t{i} = tgoal';

		% Remove reflectivity if radar was scanning  (elv <= 89 deg)
		Zcell{i}(:,elv{i}<=89) = nan;
        % Remove other measurements if radar was scanning
        VELcell{i}(:,elv{i}<=89) = nan;
        RMScell{i}(:,elv{i}<=89) = nan;
        LDRcell{i}(:,elv{i}<=89) = nan;

        % Convert status nan to 0
        status{i}(isnan(status{i})) = 0;
		% Add status flag for scanning times
		status{i}(elv{i}<=89) = 3;

        % Remove signals from cloud beards (Z < -50 dBZ)
		VELcell{i}(Zcell{i}<dbz_threshold) = nan;
		RMScell{i}(Zcell{i}<dbz_threshold) = nan;
		LDRcell{i}(Zcell{i}<dbz_threshold) = nan;
        % Remove cloud beard signals (< -50 dBZ)
		Zcell{i}(Zcell{i}<dbz_threshold) = nan;

    % If no radar files exist and this is the first loop iteration
	elseif i==1
		% Display
        disp(['File not found: ' datestr(dayvector(i), 'yymmdd')])

        % Create emtpy arrays for:
        % reflectivity
		Zcell{i} = nan(1,1);
        % other
		VELcell{i} = nan(1,1);
		RMScell{i} = nan(1,1);
		LDRcell{i} = nan(1,1);
        % time
		t{i} = (dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59))';
        % range
		h{i} = nan(1,1);
        % status
		status{i} = zeros(1, length(t{i})) + 2; % Set to 2 if file not found
        % wind
		wind{i} = zeros(1, length(t{i}));

        % 		redofirst = true;

        % Take note that dimensions have to be redone later
		redodims(i) = 1;

    % If no radar files exist and this is the not first loop iteration
    else
        % Display
		disp(['File not found: ' datestr(dayvector(i), 'yymmdd')])

        % If the range variable in the previous loop has been read and is
        % longer than one
		if length(h{i-1})>1
            % Generaty empty reflectivity array with range resolution from
            % previous day and 6*60*24 measurements per day (every ten
            % seconds)
			Zcell{i} = nan(length(h{i-1}), 60*6*24);
            % other Variables
			VELcell{i} = nan(length(h{i-1}), 60*6*24);
			RMScell{i} = nan(length(h{i-1}), 60*6*24);
			LDRcell{i} = nan(length(h{i-1}), 60*6*24);

        % If the range variable in the previous loop is only one entry
        else
            % Set one nan value to reflectivity variable
			Zcell{i} = nan(1,1);
            % other Variables
			VELcell{i} = nan(1,1);
			RMScell{i} = nan(1,1);
			LDRcell{i} = nan(1,1);
            % Take note that dimensions have to be redone later
			redodims(i) = 1;
        end

        % Create time array for this day with an entry every 10 seconds
		t{i} = (dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59))';
        % Copy range values from previous day
		h{i} = h{i-1};
        % Set status to 2 to indicate that no data file has been found
		status{i} = zeros(1, length(t{i})) + 2; % Set to 2 if file not found
        % Create empty array in wind variable
		wind{i} = zeros(1, length(t{i}));
	end

end % for i=1:length(dayvector)

% Get indices of dates where dimensions have to be redone
redodims = find(redodims);

% If dimensions have to be redone
if ~isempty(redodims)

    % Get size of each cell of reflectivity
	size_Z = cell2mat(cellfun(@size, Zcell, 'uni', false));

    % Find first cell with matrix and not only one value, i.e. first cell
    % with "real" data
	ind_firstNonNan = find(size_Z(:,1) ~= 1, 1, 'first');

    % Loop all days where dimensions have to be redone
	for i=1:length(redodims)

        % Create array of nans of the size of the first non-nan cell
		Zcell{redodims(i)} = nan(size(Zcell{ind_firstNonNan}));
        % other variables
        VELcell{redodims(i)} = nan(size(Zcell{ind_firstNonNan}));
        RMScell{redodims(i)} = nan(size(Zcell{ind_firstNonNan}));
        LDRcell{redodims(i)} = nan(size(Zcell{ind_firstNonNan}));
        % Create vector of nans of the length of the first non-nan cell
		h{i} = h{ind_firstNonNan};
	end
end

%% Adjust height dimensions

% Indices of height above set limit
ind_aboveLimit = cellfun(@(x) x>height_limit, h, 'uni', false);
% Only keep values below height limit
Zcell = cellfun(@(x, y) x(~y, :), Zcell, ind_aboveLimit, 'uni', false);
% other variables
VELcell = cellfun(@(x, y) x(~y, :), VELcell, ind_aboveLimit, 'uni', false);
RMScell = cellfun(@(x, y) x(~y, :), RMScell, ind_aboveLimit, 'uni', false);
LDRcell = cellfun(@(x, y) x(~y, :), LDRcell, ind_aboveLimit, 'uni', false);
% Only keep height below limit
h = cellfun(@(x, y) x(~y), h, ind_aboveLimit, 'uni', false);

%%%%% Make sure that all range vectors are the same %%%%%%
% Get size of range entries in cell
s = cell2mat(cellfun(@size, h, 'uni', false));

% Get unique lengths of range vector
unique_lengths = unique(s(:,1));

% Loop unique lengths of range vector
for i=1:length(unique_lengths)

    % Get index of all entries with current range vector length
    ind = s(:,1)==unique_lengths(i);

    % Convert range cells with current range vector to matrix
    height_matrix = cell2mat(h(ind)');

    % Calculate differences between entries
    differences = sum(diff(height_matrix,1,2));

    % Check if differences are zero, i.e. all range vectors are the same
    if sum(differences)~=0
%         index_diff = find(differences~=0);
		error('Height dimensions don''t match')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert data

% Get length of height arrays
h_length = cell2mat(cellfun(@(x) size(x, 1), h, 'uni', false));
% Get number of unique height arrays
uniq_h = unique(h_length);
% If more than one height array exists, find the dominant one and disregard the rest
if length(uniq_h)~=1
	for j=1:length(uniq_h)
		sum_h_length(j) = sum(h_length==uniq_h(j));
	end
	ind_discard = h_length==uniq_h(sum_h_length==min(sum_h_length));
	ind_keep = h_length==uniq_h(sum_h_length==max(sum_h_length));

	h{ind_discard} = h{ind_keep(1)};
	% Generaty empty reflectivity array with range resolution from
	% previous day and 6*60*24 measurements per day (every ten
	% seconds)
	Zcell{ind_discard} = nan(length(h{ind_keep(1)}), 60*6*24);
	% other Variables
	VELcell{ind_discard} = nan(length(h{ind_keep(1)}), 60*6*24);
	RMScell{ind_discard} = nan(length(h{ind_keep(1)}), 60*6*24);
	LDRcell{ind_discard} = nan(length(h{ind_keep(1)}), 60*6*24);

	% Create time array for this day with an entry every 10 seconds
	t{ind_discard} = (dayvector(ind_discard):1/24/60/6:dayvector(ind_discard)+datenum(0, 0, 0, 23, 59, 59))';

	% Set status to 2 to indicate that no data file has been found
	status{ind_discard} = zeros(1, length(t{ind_discard})) + 2; % Set to 2 if file not found
	% Create empty array in wind variable
	wind{ind_discard} = zeros(1, length(t{ind_discard}));
end

% Convert cell of ranges to matrix
hmat = cell2mat(h');
% Calculate differences in range gates from day to day
hdiff = diff(hmat, 1, 2);

% Error if differences between range gates vary by more than 1 cm per day
if sum(sum(abs(hdiff)>0.01))>0
	error('Something''s wrong with the range gates, please check!')
else % use first entry for all
	height = h{1};
end

% Concatenate time
time = cell2mat(t);
% Convert back to unix time
time = sdn2unixtime(time);

% Concatenate reflectivity
Z = cell2mat(Zcell');
% other variables
VEL = cell2mat(VELcell');
RMS = cell2mat(RMScell');
LDR = cell2mat(LDRcell');

% Concatenate status
status = cell2mat(status')';

% Concatenate wind
wind = cell2mat(wind')';
wind = reshape(wind, [], 1);

% ind_windnan = isnan(wind);

% Find nan entries in wind data
wind_missing = zeros(length(wind),1);
% Loop all nan entries (not pretty but hopefully working)
for i=1:length(wind)

    % If wind entry is nan
	if isnan(wind(i))

        % If current loop iteration is after first day, i.e. 6*60*24 steps
        % after the first
		if i - (6*60*24)>0

            % Take wind speed as average of previous 24 hour period
			wind(i) = mean(wind(i - (6*60*24):i-1), 'omitnan');
            % Set flag that wind was missing
			wind_missing(i) = 1;

        % If current loop iteration is within first day
        else

            % Set wind speed to 8 m/s
			wind(i) = 8;
            % Set flag that wind was missing
			wind_missing(i) = 1;
		end
	end
end

%% Save data

% Display
disp('Saving data')

if ~exist(outpath)
    mkdir(outpath)
end
% Save the data to temporary mat file
save([outpath 'Z_' radarname '_' radarrange '_' start_date '-' end_date '.mat'],...
			'Z', 'LDR', 'VEL', 'RMS', 'date','Zcell', 'RMScell', 'VELcell', 'LDRcell',...
            'height','radarname', 'radarrange','time','t','h','status', 'wind',...
            'wind_missing', '-v7.3')

% Display
disp('Concatenated data saved')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function dataout = fillData(timevec_out, timevec_in, data)

    % Get index of times to unify times
	ind_time = get_indTime(timevec_out, timevec_in);

    % Find the dimension that is not the time
	notimedim = find(size(data)~=length(timevec_in));

    % Create output data array with desired dimensions
	dataout = nan(size(data, notimedim), length(timevec_out));

    % If data is one-dimensional
	if sum(size(data)==1)
        % Put data into output array according to time indices
	    dataout(:, ~isnan(ind_time)) = data(ind_time(~isnan(ind_time)));

    % If data is two-dimensional
    else
        % Put data into output array according to time indices
	    dataout(:, ~isnan(ind_time)) = data(:, ind_time(~isnan(ind_time)));
	end
end

end

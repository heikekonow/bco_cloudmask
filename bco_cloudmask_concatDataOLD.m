function bco_cloudmask_concatData(start_date, end_date)

a = 0; b = 0; c = 0; d = 0; f = 0; g = 0; h = 0; j = 0; k = 0; l = 0; m = 0; n = 0; o = 0; p = 0; q = 0; r = 0; s = 0; t = 0; u = 0; v = 0; w = 0; x = 0; y = 0; z = 0;

height_limit = 13500;

filepath = '/pool/OBS/ACPC/MBR2/Level_1/';
outpath = '/scratch/local1/m300512/bco_concat/';
tmppath = '/scratch/local1/m300512/data_tmp/';
windpath = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/I_Meteorology_2m/';

year = start_date(1:4);

disp(['Starting processing for ' year])
disp('Concatenating')

folders = listFiles([filepath year '*'],'full');
wind_folders = listFiles([windpath year '*'],'full');

for f=1:length(folders)
	filepaths{f} = listFiles([folders{f} '/' '*.nc*'],'full');
end
for f=1:length(wind_folders)
	wind_filepaths{f} = listFiles([wind_folders{f} '/' '*.nc*'],'full');
end

files = vertcat(filepaths{:});
wind_files = vertcat(wind_filepaths{:});


% dayvector = datenum(str2num(year), 1, 1):datenum(str2num(year), 12, 31);
dayvector = datenum(start_date, 'yyyymmdd'):datenum(end_date, 'yyyymmdd');

Zcell = cell(length(dayvector),1);
t = cell(length(dayvector),1);
h = cell(length(dayvector),1);
date = cell(length(dayvector),1);
status = cell(length(dayvector),1);

redodims = zeros(1, length(dayvector));

%% Read data

for i=1:length(dayvector)
	foundfiles = strfind(files, datestr(dayvector(i), 'yymmdd'));
	ind_foundfiles = cellfun(@(x) ~isempty(x), foundfiles);

	wind_foundfiles = strfind(wind_files, datestr(dayvector(i), 'yymmdd'));
	ind_wind_foundfiles = cellfun(@(x) ~isempty(x), wind_foundfiles);

	% Make time vector with 10 second intervals (needed first time in if clause)
	tgoal = dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59);

	if [sum(ind_foundfiles)==1]

		filename = files{ind_foundfiles};
		disp(filename)

		lookforunzipped = listFiles([tmppath filename(36:end-4)]);
		% Check if files are zipped and unzip
        if isempty(lookforunzipped)
			disp('unzipping')
			eval(['! bunzip2 -c ' filename ' > ' tmppath filename(36:end-4)])
			% Rename file in list (remove extension .bz2)
			files{ind_foundfiles} = [tmppath filename(36:end-4)];
			% Check if file can be found
			if isempty(listFiles(files{ind_foundfiles}))
				error(['Error: Something''s gone wrong with renaming of unzipped file: ' files{i}])
			end
		else
			files{ind_foundfiles} = [tmppath filename(36:end-4)];
        end

        % same for wind data
		if [sum(ind_wind_foundfiles)==1]
			wind_filename = wind_files{ind_wind_foundfiles};
	        lookforunzipped = strcmp(wind_filename(end-1:end), 'nc');
	        if ~lookforunzipped & isempty(listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']))
				disp('unzipping wind')
				if ~isempty(strfind(wind_filename, 'bz'))
					eval(['! bunzip2 -c ' wind_filename ' > ' tmppath wind_filename(70:end-4)])
					% Rename file in list (remove extension .bz2)
					wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end-4)];
				else
					eval([' ! cp ' wind_filename ' ' tmppath])

					% Rename file in list (set new path)
					wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end)];
				end

				% Check if file can be found
				if isempty(listFiles(wind_files{ind_foundfiles}))
					error(['Error: Something''s gone wrong with renaming of unzipped file: ' files{i}])
				end
			elseif isempty(listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']))
				disp('copying wind')

				eval([' ! cp ' wind_filename ' ' tmppath])

				% Rename file in list (set new path)
				wind_files{ind_wind_foundfiles} = [tmppath wind_filename(70:end)];

			else
				% Rename file in list (set new path)
				wind_out = listFiles([tmppath 'Meteorology*' datestr(dayvector(i), 'yymmdd') '*']);
				wind_files{ind_wind_foundfiles} = [tmppath wind_out{1}];
	        end

			% Read wind data if the file exists
	        windread = ncread(wind_files{ind_wind_foundfiles}, 'VEL');
	        twindread = unixtime2sdn(ncread(wind_files{ind_wind_foundfiles}, 'time'));
		else
			windread = nan(length(tgoal),1);
			twindread = tgoal;
		end



		% Read data
		if ncVarInFile(files{ind_foundfiles},'Zf')
			Zread= ncread(files{ind_foundfiles},'Zf');
		elseif ncVarInFile(files{ind_foundfiles},'Z')
			Zread= ncread(files{ind_foundfiles},'Z');
		else
			error(['Variable Z or Zf not found in file' files{ind_foundfiles}])
		end
		tread = unixtime2sdn(ncread(files{ind_foundfiles},'time'));
        statusread = ncread(files{ind_foundfiles},'status');


        Zcell{i} = fillData(tgoal, tread, Zread);
        t{i} = tgoal';

		h{i} = ncread(files{ind_foundfiles},'range');
		status{i} = fillData(tgoal, tread, statusread);
        status{i}(isnan(status{i})) = 0;
        wind{i} = fillData(tgoal, twindread, windread);

		% Remove cloud beard signals (< -50 dBZ)
		Zcell{i}(Zcell{i}<-50) = nan;
	elseif i==1
		disp(['File not found: ' datestr(dayvector(i), 'yymmdd')])
		Zcell{i} = nan(1,1);
		t{i} = (dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59))';
		h{i} = nan(1,1);
		status{i} = zeros(1, length(t{i})) + 2; % Set to 2 if data filled
		wind{i} = zeros(1, length(t{i}));
		redofirst = true;
		redodims(i) = 1;
	else
		disp(['File not found: ' datestr(dayvector(i), 'yymmdd')])
		if length(h{i-1})>1
			Zcell{i} = nan(length(h{i-1}), 60*6*24);
		else
			Zcell{i} = nan(1,1);
			redodims(i) = 1;
		end
		t{i} = (dayvector(i):1/24/60/6:dayvector(i)+datenum(0, 0, 0, 23, 59, 59))';
		h{i} = h{i-1};
		status{i} = zeros(1, length(t{i})) + 2; % Set to 2 if data filled
		wind{i} = zeros(1, length(t{i}));
	end

%		imagesc(t{i},h{i},Zcell{i}); set(gca,'YDir','normal')
%		xdate
%		date{i} = files{i}(end-8:end-3);
%		title(date{i})
end

redodims = find(redodims);

if ~isempty(redodims)
	size_Z = cell2mat(cellfun(@size, Zcell, 'uni', false));
	ind_firstNonNan = find(size_Z(:,1) ~= 1, 1, 'first');
	for i=1:length(redodims)
		Zcell{redodims(i)} = nan(size(Zcell{ind_firstNonNan}));
		h{i} = h{ind_firstNonNan};
	end
end

%% Adjust height dimensions
% Indices of height above set limit
ind_aboveLimit = cellfun(@(x) x>height_limit, h, 'uni', false);
% Only keep values below height limit
Zcell = cellfun(@(x, y) x(~y, :), Zcell, ind_aboveLimit, 'uni', false);
% Only keep height below limit
h = cellfun(@(x, y) x(~y), h, ind_aboveLimit, 'uni', false);

%%%%%%%%%%%
s = cell2mat(cellfun(@size, h, 'uni', false));

unique_lengths = unique(s(:,1));

for i=1:length(unique_lengths)
    ind = s(:,1)==unique_lengths(i);

    height_matrix = cell2mat(h(ind)');

    differences = sum(diff(height_matrix,1,2));

    if sum(differences)~=0
        index_diff = find(differences~=0);
		error('Height dimensions don''t match')
	end
end
%%%%%%%%%%%

%% Convert data
% Convert cell of ranges to matrix
hmat = cell2mat(h');
% Calculate differences in range gates from day to day
hdiff = diff(hmat,2);

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

% Concatenate status
status = cell2mat(status')';

% Concatenate wind
wind = cell2mat(wind')';
wind = reshape(wind, [], 1);

% Find nan entries in wind data
ind_windnan = isnan(wind);
wind_missing = zeros(length(wind),1);
% Loop all nan entries (not pretty but hopefully working)
for i=1:length(wind)
	if isnan(wind(i))
		wind(i) = mean(wind(i - (6*60*24):i-1), 'omitnan');
		wind_missing(i) = 1;
	end
end

%% Save data

disp('Saving data')

% save([outpath 'Z_' year '-01.mat'],'Z','date','Zcell','height','time','t','h','status', 'wind','-v7.3')
save([outpath 'Z_' start_date '-' end_date '.mat'],'Z','date','Zcell','height',...
			'time','t','h','status', 'wind', 'wind_missing', '-v7.3')

disp('Concatenated data saved')

%% Functions
function dataout = fillData(timevec_out, timevec_in, data)

	ind_time = get_indTime(timevec_out, timevec_in);

	notimedim = find(size(data)~=length(timevec_in));

	% dataout = nan(min(size(data)), length(timevec_out));
	dataout = nan(size(data, notimedim), length(timevec_out));

	if sum(size(data)==1)
	    dataout(:, ~isnan(ind_time)) = data(ind_time(~isnan(ind_time)));
	else
	    dataout(:, ~isnan(ind_time)) = data(:, ind_time(~isnan(ind_time)));
	end
end

end

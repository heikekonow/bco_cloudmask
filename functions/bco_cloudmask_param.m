% 	Code for calculating cloud parameters from segmented cloud objects.
% 	Input variables:
%		- start_date	string with first date to process (yyyymmdd format)
%		- end_date		string with last date to process (yyyymmdd format)
%		- radarname		string with radar name (MBR or KATRIN)
%		- radarrange	string with radar height range (e.g. '155m-18m')
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020



function bco_cloudmask_param(start_date, end_date, radarname, radarrange)

% Necessary, since this code has been copied from flight data processing where
% multiple flights were looped over. This should be removed in the future.
flightdate = 1;
% Set to true if test figures should be produced
figures = false;

% Inform user
disp('Calculating parameters')

%% Load data %%%%%%%%%%%%%%%
% Set path to data to be read
filepath = ['/scratch/local1/m300512/bco_concat/Z_' radarname '_' radarrange '_' start_date '-' end_date '_closed.mat'];

% Read data
datastruct = load(filepath);

% Rename variables
conComp = datastruct.con_comp;
time_t = datastruct.time;
height_t = datastruct.height;
wind = datastruct.wind;
con_comp{1} = conComp;
time{1} = time_t;
height{1} = height_t;
clear conComp time_t height_t


%% Calculate parameters %%%%%%%%%%%%%%%

% Preallocate arrays
windSpeed = cell(length(flightdate),1);
maskIndClouds = cell(length(flightdate),1);
averageWindSpeed = cell(length(flightdate),1);
cloudLength = cell(length(flightdate),1);
cloudDepth = cell(length(flightdate),1);
cloudBase = cell(length(flightdate),1);
cloudTop = cell(length(flightdate),1);
cloudStartTime = cell(length(flightdate),1);
cloudEndTime = cell(length(flightdate),1);

% Loop flights (unnecessary for BCO data)
for i=1:length(flightdate)

    % Number of clouds
    numberOfClouds = con_comp{i}.NumObjects;

    % Generate matrices of height and time values
    height_mat = repmat(height{i},1,length(time{i}));
    time_mat = repmat(time{i}', length(height{i}),1);

    % Calculate wind speed in all heights using power law
    wind_2m_mat = repmat(wind', size(height_mat,1), 1);
    windSpeed{i} = wind_2m_mat .* (height_mat ./ 2).^0.11;

    % Loop individual clouds
    for j=1:numberOfClouds

        % Lowest cloud height
        minHeight = min(height_mat(con_comp{i}.PixelIdxList{j}));
        % Heighest cloud height
        maxHeight = max(height_mat(con_comp{i}.PixelIdxList{j}));

        % First cloud time step
        minTime = min(time_mat(con_comp{i}.PixelIdxList{j}));
        % Last cloud time step
        maxTime = max(time_mat(con_comp{i}.PixelIdxList{j}));

        if figures
            if j==1
                figure; set(gcf,'Position',[-875 486 828 595])
                cm = brewermap(10000,'Set1');
                colormap(cm)
                addWhiteToColormap
            end
            imagesc(time{i},height{i},maskIndClouds{i}{j})
            set(gca,'XLim',[minTime-10 maxTime+10],'YLim',[minHeight-70 maxHeight+70],...
                    'YDir','normal')
            grid on
            pause;
        end

        % Convert linear cloud indices to rows/columns
        [~,col] = ind2sub([length(height{i}),length(time{i})],con_comp{i}.PixelIdxList{j});

        %%% Cloud length %%%
        % Average ground speed for this cloud
        averageWindSpeed{i}(j) = mean(windSpeed{i}(con_comp{i}.PixelIdxList{j}));
        % Cloud length in time space; + 1 because a one second long cloud
        % fills one profile, so this is then 200 m long
        cloudLengthTime = maxTime-minTime + 1;                          % seconds
        % Cloud length from time and velocity: l = v * t
        cloudLength{i}(j) = averageWindSpeed{i}(j) .* cloudLengthTime;  % meters

        %%% Cloud depth %%%
        cloudDepth{i}(j) = maxHeight - minHeight;

        %%% Cloud base %%%
        cloudBase{i}(j) = minHeight;

        %%% Cloud top %%%
        cloudTop{i}(j) = maxHeight;

        %%% Cloud time steps
        cloudStartTime{i}(j) = minTime;
        cloudEndTime{i}(j) = maxTime;
     end

end

%% Save data %%%%%%%%%%%%%%%

Z = datastruct.Z;
LDR = datastruct.LDR;
VEL = datastruct.VEL;
RMS = datastruct.RMS;
date = datastruct.date;
numMask = datastruct.numMask;
size_se = datastruct.size_se;
se = datastruct.se;
type_se = datastruct.type_se;
status = datastruct.status;
wind_missing = datastruct.wind_missing;
cloudMask = datastruct.cloudMask;
cloudMask_closed = datastruct.cloudMask_closed;

save([filepath(1:end-4) '_concomp.mat'],'Z', 'LDR', 'VEL', 'RMS','date','height','time*',...
    'averageWindSpeed','cloud*','con_comp','windSpeed','height*',...%'maskIndClouds',...
    'numMask','numberOfClouds','size_se','se','type_se','status', 'wind_missing',...
    '-v7.3')

disp('Parameters saved')

function bco_cloudmask_param(start_date, end_date, radarname, radarrange)

% clear; % close all

flightdate = 1;
figures = 0;

disp('Calculating parameters')

%% Load data
filepath = ['/scratch/local1/m300512/bco_concat/Z_' radarname '_' radarrange '_' start_date '-' end_date '_closed.mat'];

load(filepath)
% Rename variable
filepath = ['/scratch/local1/m300512/bco_concat/Z_' radarname '_' radarrange '_' start_date '-' end_date '_closed.mat'];

conComp = con_comp;
time_t = time;
height_t = height;
clear con_comp time height
con_comp{1} = conComp;
time{1} = time_t;
height{1} = height_t;
clear conComp time_t height_t


%% Calculate parameters

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

% Loop flights
for i=1:length(flightdate)


%    groundSpeed{i} = repmat(8,1,length(time{i}));
%     windSpeed{i} = wind;


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
        % Generate masks for each individual cloud
%        maskIndClouds{i}{j} = zeros(length(height{i}),length(time{i}));
%        maskIndClouds{i}{j}(con_comp{i}.PixelIdxList{j}) = 1;

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
%        [~,col] = ind2sub(size(maskIndClouds{i}{j}),con_comp{i}.PixelIdxList{j});
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

%% Save data

save([filepath(1:end-4) '_concomp.mat'],'Z','date','height','time*',...
    'averageWindSpeed','cloud*','con_comp','windSpeed','height*',...%'maskIndClouds',...
    'numMask','numberOfClouds','size_se','se','type_se','status', 'wind_missing',...
    '-v7.3')
%save([filepath(1:end-4) '_concomp.mat'],'Z','date','Zcell','height','time*','t','h',...
%    'averageGroundSpeed','cloud*','con_comp','groundSpeed','height*',...%'maskIndClouds',...
%    'numMask','numberOfClouds','size_se','se','type_se','status',...
%    '-v7.3')

disp('Parameters saved')

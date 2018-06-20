function bco_cloudmask_mask(radarname, radarrange, start_date, end_date)
% function bco_cloudmask_mask(start_date, end_date)
    % radarname{i}, unique_height{i}{j}, dates{i,j}(ind_years,:)

% start_date = dates(1,:);
% end_date = dates(end,:);
%% Load data

% clear; %close all

% filepath = ['/scratch/local1/m300512/bco_concat/Z_' start_date '-' end_date '.mat'];
filepath = ['/scratch/local1/m300512/bco_concat/Z_' radarname '_' radarrange '_' start_date '-' end_date '.mat'];

disp('Generating cloud mask')

load(filepath, 'Z')

%% Generate cloud mask
cloudMask = zeros(size(Z));

cloudMask(~isnan(Z)) = 1;

%% Image closing

disp('Closing image')

% Set type and size of structuring element
type_se = 'rectangle';
size_se = [5 2];

% Define structuring Element
se = strel(type_se,size_se);

% Close cloud mask
cloudMask_closed = imclose(cloudMask,se);

%% Connected components labeling from initial radar reflectivity, now after closing

disp('Connected components labeling')
% Apply connected components labeling (and remove clouds that are
% smaller than 4 pixels)
[con_comp,numMask] = radar_connectedClouds(cloudMask_closed);



%% Save data

disp('Saving data')

save([filepath(1:end-4) '_closed.mat'], '-v7.3')

clear Z cloudMask_closed con_comp numMask

load(filepath, 'height', 'time', 'status', 'wind', 'date', 'wind_missing', 'LDR', 'VEL', 'RMS')

save([filepath(1:end-4) '_closed.mat'], '-append')

disp('Cloudmask data saved')

% 	Code that generates cloud masks from radar reflectivity data.
%
% 	Input variables:
%		- radarname		string with radar name (MBR or KATRIN)
%		- radarrange	string with radar height range (e.g. '155m-18m')
%		- start_date	string with first date to process (yyyymmdd format)
%		- end_date		string with last date to process (yyyymmdd format)
%       - minSize       minimum size for objects, everything smaller will
%                       be removed as clutter
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020

function bco_cloudmask_mask(radarname, radarrange, start_date, end_date, minSize, tmppath)


%% Load data %%%%%%%%%%%%%%%%%%%

% Specify path to temporary files
filepath = [tmppath 'Z_' radarname '_' radarrange '_' start_date '-' end_date '.mat'];

% Inform user
disp('Generating cloud mask')

% Read data from temporary file
load(filepath, 'Z')

%% Generate cloud mask (0 = no cloud, 1 = cloud)
cloudMask = zeros(size(Z));
cloudMask(~isnan(Z)) = 1;


%% Image closing %%%%%%%%%%%%%%%%%%%

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
[con_comp,numMask] = radar_connectedClouds(cloudMask_closed, minSize);



%% Save data %%%%%%%%%%%%%%%%%%%

disp('Saving data')

save([filepath(1:end-4) '_closed.mat'], '-v7.3')

clear Z cloudMask_closed con_comp numMask

load(filepath, 'height', 'time', 'status', 'wind', 'date', 'wind_missing', 'LDR', 'VEL', 'RMS')

save([filepath(1:end-4) '_closed.mat'], '-append')

disp('Cloudmask data saved')

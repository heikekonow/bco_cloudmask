function [con_comp,numMask] = radar_connectedClouds(cloudMask)

% Preallocate
numMask = nan(size(cloudMask));

% Replace nan with zeros
cloudMask(isnan(cloudMask)) = 0;

% Get structure with connected components informations:
% cc.Connectivity - Connectivity of the connected components (objects)
% cc.ImageSize - Size of cloudMask
% cc.NumObjects - Number of connected components
% cc.PixelIdxList - 1-by-NumObjects cell array where the k-th element in the 
%                   cell array is a vector containing the linear indices of 
%                   the pixels in the k-th object. 
con_comp = bwconncomp(cloudMask);

%% Filter out clouds that are smaller than 4 pixels
% Get size of single clouds
comp_size = cell2mat(cellfun(@length,con_comp.PixelIdxList,'UniformOutput',false));
% Define index for clouds smaller than 4 pixels
ind_toosmall = comp_size<4;
% Delete clouds that are too small
con_comp.PixelIdxList(ind_toosmall) = [];
% Remove number of too small clouds
con_comp.NumObjects = con_comp.NumObjects-sum(ind_toosmall);

%% Add cloud id
% Loop all found Objects
for i=1:con_comp.NumObjects
    % assign each found cloud a consecutive number
    numMask(con_comp.PixelIdxList{i}) = i;
end
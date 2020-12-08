%% [con_comp,numMask] = radar_connectedClouds(cloudMask)
% 	Function to perform connected components labelling on radar cloud
% 	masks.
%
% 	Input variables:
%		- cloudMask 	binary cloud mask (0 = no cloud, 1 = cloud)
%
%   Output variables:
%       - con_comp      structure connected components informations
%           - Connectivity  Connectivity of the connected components (objects)
%           - ImageSize     Size of cloudMask
%           - NumObjects    Number of connected components
%           - PixelIdxList 	1-by-NumObjects cell array where the k-th element in the 
%                           cell array is a vector containing the linear indices of 
%                           the pixels in the k-th object. 
%       - numMask       array of the same size as cloudMask where each
%                       cloudy pixel is identified by their object ID
%
%	contact: Heike Konow, heike.konow@uni-hamburg.de
%	last revision: Dec 2020

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
function configstruct = config

% This is the config file for BCO cloud mask processing
%
% 1. change the name of this file to config.m
% 
% 2. Set the paths to data directory below

% Set path to radar data
configstruct.path = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/Ka-Band/MBR2/10s/';
% Set path to output files
configstruct.outpath = '/pool/OBS/ACPC/MBR2/cloudmask/bco_object_cloudmask/';
% Set path for temporary files
configstruct.outpathtmp = '/scratch/local1/m300512/bco_concat/';
% Set path for temporary files
configstruct.tmppath = '/scratch/local1/m300512/data_tmp/';
% Set path to wind data
configstruct.windpath = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/I_Meteorology_2m/';
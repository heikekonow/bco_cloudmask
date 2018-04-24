% Clean up
clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set date range to work with
% start_date_all = '20150501';
% end_date_all = '20180331';
% start_date_all = '20160101';
start_date_all = '20150601';    % Achtung, MBR erst ab Juni 2015 nehmen, davor komische Hoehen
% end_date_all = '20160431';
end_date_all = '20150610';
% end_date_all = '20160530';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare dates %%%%%%%%%%%%%%%%%%%%%%%%
% Make vector with all days between start_date and end_date
dayvector = datenum(start_date_all, 'yyyymmdd'):datenum(end_date_all, 'yyyymmdd');
% Convert to datevector
dayvector_vec = datevec(dayvector);
% Only year and months
yyyymm = str2num([num2str(dayvector_vec(:,1)) num2str(dayvector_vec(:,2), '%02d')]);
% yyyymm = str2num([num2str(dayvector_vec(:,1))]);
% Year, months, and days
yyyymmdd = [num2str(dayvector_vec(:,1)) num2str(dayvector_vec(:,2), '%02d') num2str(dayvector_vec(:,3), '%02d')];
% Look for unique year/month combinations
months_use = unique(yyyymm);

% Loop unique months
for i=1:length(months_use)

    disp(['>>>>>>>> Processing ' num2str(months_use(i)) ' <<<<<<<<'])
    % Find indices in entire date vector for current month
    ind_months = strncmp(cellstr(num2str(months_use(i))), cellstr(yyyymmdd), 4);

    % Make vector with days to use for this
    days_use = yyyymmdd(ind_months,:);

    % Get first and last day to use in functions
    start_date = days_use(1,:);
    end_date = days_use(end,:);

    %% Actual processing %%%%%%%%%%%%%%%%%%%%%%%%

    % Concatenate data
    bco_cloudmask_concatData(start_date, end_date)

    % Generate cloud mask
    bco_cloudmask_mask(start_date, end_date)

    % Caclulate cloud parameter
    bco_cloudmask_param(start_date, end_date)

    % Save data to netcdf
    bco_cloudmask_save2netcdf(start_date, end_date)
end

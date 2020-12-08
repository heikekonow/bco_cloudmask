%% get_indTime(uniTime,instrTime, varargin)
%   get_indTime - Get indices of instrument time to transfer to uniform time
%
%   Syntax:  indTime = get_indTime(uniTime,instrTime, varargin)
%
%   Input variables:
%       - uniTime   uniform time to which instrument time should be transfered
%       - instrTime instrument measurement time
%       - varargin  instrument name (only set to be relevant if instrument
%                   is 'bahamas'
%
%   Output variables:
%       - indTime   index array for instrument time to transfer to uniform time
%
%   contact: Heike Konow, heike.konow@uni-hamburg.de
%   March 2016; Last revision: June 2016

%%
function indTime = get_indTime(uniTime,instrTime, varargin)

%------------- BEGIN CODE --------------

% Value for one second
oneSecond = 1/24/60/60;

% If no extra information is provided
if nargin==2
    % Preallocate
    indTime = nan(1,length(uniTime));
    % Loop all uniform time steps
    for i=1:length(uniTime)
        % Calculate difference between each time step and all others
        absDifference = abs(instrTime-uniTime(i));

        % Find index of minimal time difference, i.e. the time step closest
        % to selected uniform time step
        indMinim(i) = find(absDifference==min(absDifference),1,'first');

        % If this is not the first time step and index is not the same as
        % the one before
        if i>1 && (indMinim(i)~=indMinim(i-1))% && ~isnan(indTime(i-1)))
            % Save time index
            indTime(i) = indMinim(i);
        % If this is the first time step
        elseif i==1
            % Save time index
            indTime(i) = indMinim(i);
        end
    end

% If additional info is given
elseif nargin==3
    % Check Name
    if strcmp(varargin{1},'bahamas')
        % Find first and last time step in Bahamas times that match uniform
        % time
        a = find(uniTime<instrTime(1),1,'last');
        b = find(uniTime>instrTime(end),1,'first');
        
        % Get length of relevant time interval
        t = b-a+1;
        
        % Preallocate
        indTimeUni = nan(t,1);
        indTimeInstr = nan(t,1);
        % Initialize counter
        k = 1;
        
        % Loop time
        for i=a:b
            
            % Find minimal time deviation for current time step between
            % uniform time and Bahamas time
            absDifference = abs(instrTime-uniTime(i));
            indMinim(k) = find(absDifference==min(absDifference),1,'first');

            % If this is not the first time step
            if k>1 && (indMinim(k)~=indMinim(k-1)) && absDifference(indMinim(k))<oneSecond% && ~isnan(indTime(i-1)))
                indTimeInstr(k) = indMinim(k);
                indTimeUni(k) = i;
                
            % If this is the first time step
            elseif k==1
                indTimeInstr(k) = indMinim(k);
                indTimeUni(k) = i;
            end
            % Increase counter
            k = k+1;
        end
    else
        error('Only Bahamas is implemented. Use argument ''bahamas''.')
    end
else
     error('Wrong number of input arguments')
end


end
%------------- END OF CODE --------------

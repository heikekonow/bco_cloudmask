%% get_indTime
%   get_indTime - One line description of what the function or script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Example:
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also:
%
%   Author: Dr. Heike Konow
%   Meteorological Institute, Hamburg University
%   email address: heike.konow@uni-hamburg.de
%   Website: http://www.mi.uni-hamburg.de/
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
        a = find(uniTime<instrTime(1),1,'last');
        b = find(uniTime>instrTime(end),1,'first');
        % length of relevant time interval
        t = b-a+1;

        indTimeUni = nan(t,1);
        indTimeInstr = nan(t,1);
        k = 1;
        for i=a:b
            absDifference = abs(instrTime-uniTime(i));

            indMinim(k) = find(absDifference==min(absDifference),1,'first');

            if k>1 && (indMinim(k)~=indMinim(k-1)) && absDifference(indMinim(k))<oneSecond% && ~isnan(indTime(i-1)))
                indTimeInstr(k) = indMinim(k);
                indTimeUni(k) = i;
            elseif k==1
                indTimeInstr(k) = indMinim(k);
                indTimeUni(k) = i;
            end
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

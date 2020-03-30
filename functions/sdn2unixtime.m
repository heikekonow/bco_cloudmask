%% sdn2unixtime
%   sdn2unixtime - Convert Matlab serial date numbers to Unix Time (seconds
%                  since 01.01.1970 00:00:00) 
%
%   Syntax:  tUnix = sdn2unixtime(tSDN)
%
%   Inputs:
%       tSDN - time as Serial Date Number
%
%   Outputs:
%       tUnix - time in seconds since 01.01.1970 00:00:00
%
%   Example: 
%       sdn2unixtime(now)
%
%       ans =
%           1434109923.0384
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: unixtime2sdn
%
%   Author: Dr. Heike Konow
%   Meteorological Institute, Hamburg University
%   email address: heike.konow@uni-hamburg.de
%   Website: http://www.mi.uni-hamburg.de/
%   January 2014; Last revision: June 2015

%%
function tUnix = sdn2unixtime(tSDN)

%------------- BEGIN CODE --------------

oneSecond = datenum(0,0,0,0,0,1);
tUnix = (tSDN - datenum(1970,1,1,0,0,0)) ./ oneSecond;

%------------- END OF CODE --------------
%% unixtime2sdn
%   unixtime2sdn - Convert Unix Time (seconds since 01.01.1970 00:00:00) to
%                   Matlab serial date numbers
%
%   Syntax:  tSDN = unixtime2sdn(tUnix)
%
%   Inputs:
%       tUnix - time in seconds since 01.01.1970 00:00:00
%
%   Outputs:
%       tSDN - time as Serial Date Number
%
%   Example: 
%       tUnix = 1434109923.0384;
%       unixtime2sdn(tUnix)
%
%       ans =
%           736127.494479611
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: sdn2unixtime
%
%   Author: Dr. Heike Konow
%   Meteorological Institute, Hamburg University
%   email address: heike.konow@uni-hamburg.de
%   Website: http://www.mi.uni-hamburg.de/
%   January 2014; Last revision: June 2015

%%
function tSDN = unixtime2sdn(tUnix)
% Convert Unix Time (seconds since 01.01.1970 00:00:00) to Matlab serial
% date number

secondsPerDay = 86400;

if isinteger(tUnix)
    tUnix = double(tUnix);
end

tSDN = tUnix ./ secondsPerDay + datenum(1970,1,1,0,0,0);
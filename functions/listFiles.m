%% listFiles
%   listFiles - listing file names in given directory
%
%   Syntax:  Files = listFiles(searchstring)
%
%   Inputs:
%       searchstring - string to search for, i.e. path; may be combined
%                      with wildcards and file extensions
%       varargin -     'fullpath' (optional) gives full path to file
%                      'last'/'latest' (optional) gives newest file version
%
%   Outputs:
%       Files - cell array with file names as strings
%
%   Example:
%       BahamasPath = '/data/share/u231/u231107/HAMP/bahamas_all/';
%       files = listFiles([BahamasPath '*.nc'])
%       files =
%             'adlr_20131210a.naswNaNs.nc'
%             'adlr_20131211a.naswNaNs.nc'
%             'adlr_20131212a.naswNaNs.nc'
%             ...
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
%   January 2014,
%   June 2015,
%   February 2017: added option for full path return
%   September 2017: added option for latest file version return
%   March 2020: remove entries '.' and '..' for folder lists
%

%%
function Files = listFiles(searchstring,varargin)

%------------- BEGIN CODE --------------

filescells = cell(size(searchstring,1),1);
for i=1:size(searchstring,1)
    % List every file that complies to searchstring in structure
    filestruct = dir(searchstring(i,:));
    % Convert structure to cell
    filestruct = struct2cell(filestruct);
    % Only the first entry (file names) is needed
    filescells{i,1} = filestruct(1,:)';
end

Files = vertcat(filescells{:});

% Remove entries for directories '.' and '..'
Files(strcmp(Files, '.')) = [];
Files(strcmp(Files, '..')) = [];

% If extra options are given
if nargin>1
    % If full path return is called
    if any(strncmp(varargin,'full',4))
        % Find slashes in search path
        ind_slashes = regexp(searchstring,'/');
        % Path is the whole path until the last slash
        path = {searchstring(1:ind_slashes(end))};
        % Match size to number of files found
        path = repmat(path,size(Files));
        % Combine path and file names
        Files = cellstr(horzcat(char(path),char(Files)));
    end

    % If latest version is called
    if any(strcmp(varargin,'last')) || any(strcmp(varargin,'latest'))

        % Preallocate
        v = zeros(length(Files), 2);
        % Loop all files
        for i=1:length(Files)
            % Analyze file names for version string
            ind_startVersion = regexp(Files{i}, '_v')+2;
            ind_fileExtension = regexp(Files{i}(ind_startVersion:end), '.nc');

            % Get short file name
            filename_short = Files{i}(ind_startVersion:ind_startVersion+ind_fileExtension-2);

            % Find index of dots
            ind_dots = regexp(filename_short, '[.]');
            % If dots were found
            if ~isempty(ind_dots)
                % Analyse file version and subversion
                versionstring = filename_short(1:ind_dots(1)-1);
                subversionstring = filename_short(ind_dots(1)+1:end);

                v(i,:) = [str2double(versionstring) str2double(subversionstring)];
            end
        end
        % Convert version to float
        v = v(:,1) + v(:,2).* .01;

        % If versions found are larger than zero
        if sum(sum(v))>0
            % Get index of highest version number
            ind_max = v==max(v);
            % Get file name of highest version number
            Files = Files{ind_max};
        else
            % Otherwise take last file
            Files = Files{end};
        end
    end
end

%------------- END OF CODE --------------

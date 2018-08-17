function resp = afpc_filemover(filename_out, dir_sec, abbrev, sys) %#ok<INUSD> legacy support

% filemover: moves files in the anova_fullpipecall
%   afpc_filemover(filename_out, dir_sec, abbrev) moves files that begin 
%   with the string in filename_out from the current directory to the one 
%   specified by abbrev. The directory that the files end up in is
%   dir_sec(1:end-1)_abbrev. If the directory does not exist, an effort is
%   made to create it. 
%
%   afpc_filemover(filename_out, dir_sec, abbrev, sys) takes an
%   additional argument that specifies the value of the error that should
%   be expected if a new directory is to be made. See function for more
%   details.
%
%   resp = sfpc_filemover(...) returns the final error message from the
%   file move. If resp == 0, then the move was successful.

% Written by Rannee Li, Jul 2015
% Update log
%  - Mar 2016: added ability to deal with spaces in path.
%   18 Jul 17:  removed use of sys, resp and all that because 'exist' is a
%               function that exists. 

% specify the name of the output directory
out_dir = [dir_sec(1:(end-1)) '_' abbrev '/'];

% deal with spaces in file path
k = strfind(out_dir, ' ');

if ~isempty(k)
    warning('Spaces detected in path name')
    for a = 1:numel(k)
        out_dir(k(a):end+1) = ['\' out_dir(k(a):end)];
    end
end

% find out if the folder exists.
if exist(out_dir, 'dir')~=7
    % directory does not exist, so make it. 
    warning('Destination directory not found. Making directory %s.', out_dir)
    system(['mkdir ' out_dir]);
end

% move the files
resp = system(['mv ' filename_out '* ' out_dir]);

if nargout < 1
    clear resp
end
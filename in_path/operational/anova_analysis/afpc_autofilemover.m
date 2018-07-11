function resp = afpc_autofilemover(filename_out, dir_sec, abbrev)

% filemover: moves files in the anova_fullpipecall
%   afpc_filemover(filename_out, dir_sec, abbrev) moves files that begin 
%   with the string in filename_out from the current directory to the one 
%   specified by abbrev. The directory that the files end up in is
%   dir_sec(1:end-1)_abbrev. If the directory does not exist, an effort is
%   made to create it. 
%
%   resp = afpc_filemover(...) returns the final error message from the
%   file move. If resp == 0, then the move was successful.

% Written by Rannee Li, Jul 2015
% Update log
%  - Mar 2016: added ability to deal with spaces in path.

% % automatically determine the expected value of sys
% sysstr = computer;
% switch sysstr
%     case 'PCWIN64'
%         warning('Windows system detected. Unsure of what the error number is.')
%     case 'GLNXA64'
%         sys = 1;
%         % in ubuntu, all errors return 1, so even if it's not a 'this
%         % directory does not exist' issue.
%     case 'MACI64'
%         sys = 64;
% end



% specify the name of the output directory
out_dir = [dir_sec(1:(end-1)) '_' abbrev '/'];

% deal with spaces in file path
k = strfind(out_dir, ' ');

for a = 1:numel(k)
    warning('Spaces detected in path name')
    out_dir(k(a):end+1) = ['\' out_dir(k(a):end)];
end


% find out if the directory exists
isthere = exist(out_dir, 'dir');

if isthere~=7
    % directory does not exist. make it.
    warning('Destination directory not found. Making directory.')
    system(['mkdir ' out_dir]);
end
    
    
resp = system(['mv ' filename_out '*.mat ' out_dir]);    


if nargout < 1
    clear resp
end
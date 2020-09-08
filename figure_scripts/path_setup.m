% Adds paths to workspace and path

%% Add paths to workspace
% directory where the data is mounted. 
% Specifpy path to your project directory.
% e.g. disk_path = 'D:/yota/'
disk_path = 'your project directory';

% set figure path to be this temporary directory that is not on the main
% script or figure paths (this is a directory where the old and broken
% scripts are kept).
fig_path = [];

% highest level directory of data
data_path = [disk_path 'data/'];

% directory with the raw files
raw_path = [data_path 'raw/'];

% path to third party toolboxes
tool_path = [disk_path 'scripts/thirdparty_toolboxes/'];
chron_dir = [tool_path 'chronux/'];

% Add path to disk
% addpath(genpath(disk_path));
%% Add repository to path
%repo_path = fullfile(mfilename('fullpath'), '../../');
[rpath, rname, rext] = fileparts(mfilename('fullpath'));
repo_path = fullfile(fpath, '../../');

% check if already on path
path_list = regexp(path, pathsep, 'split');
if ~any(ismember(repo_path, path_list))
    addpath(genpath(repo_path))
end

clear path_list
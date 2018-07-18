%% get overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% set specific variables
% find filenames
data_dir = fullfile(data_path, 'included_datasets'); 
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

%% Get data filenames
% find out if the data exists, where it is, and if it doesn't exist, make
% it.

% Since the names follow a known pattern, generate one and then generate
% the rest by repeating it.

% Find name for a single cat
cat_names = dirsinside(data_dir);
dir_sec = fullfile(data_dir, cat_names{1}, data_type, '/');
loadname = dir(fullfile(dir_sec, '*.mat'));

% Check if the files exist. If not, make them
if isempty(loadname)
    extract_fstat(data_path, data_type(1:end-8))
    loadname = dir(fullfile(dir_sec, '*.mat'));
end

% Make the struct name into a matrix
nFiles = numel(loadname);
loadname = struct2cell(loadname);
loadname = cell2mat(loadname(1, :)');
loadname = horzcat(repmat(dir_sec, [nFiles, 1]), loadname);

% duplicate for the rest of the cats
nCats = numel(cat_names);
all_names = repmat(loadname, [nCats,1]);
data_dir_ind = length(data_dir) + 2;
dir_sec_ind = length(dir_sec) + 1;
cat_length = length(cat_names{1});
for k = 2:nCats
    for m = 1:nFiles
        all_names((k-1)*nFiles+m, data_dir_ind:data_dir_ind+cat_length-1) = cat_names{k};
        all_names((k-1)*nFiles+m, dir_sec_ind:dir_sec_ind+cat_length-1) = cat_names{k};
    end
end

%% Load files and pre-process


% Calls the p_thresholder_filenames function

%% Set your selection here
% Set q here
q = 0.05;

% Run per file? false = run for all files together.
ifperfile = false;

%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Find list of files

% find a list of cats
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));
dtype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout';
fnames = cell(numel(cat_names)*2, 1);

for c = 1:numel(cat_names)
    cDir = fullfile(data_path, 'included_datasets', cat_names{c}, dtype);
    files = dir(fullfile(cDir, '*.mat'));
    for k = 1:2
        fnames{(c-1)*2 + k} = fullfile(cDir, files(k).name);
    end
end

%% Call function
p_thresholder_filenames(fnames, q, ifperfile);

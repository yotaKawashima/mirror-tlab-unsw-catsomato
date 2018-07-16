% Calls the p_thresholder function

%% Set your selection here
% Put a string with cat name. Leave empty to grab all files. 
cats = []; 

% Put an int for area. Leave empty for all areas.
areas = [];

% Set q here
q = 0.05;

% Run per file? false = run for all files together.
ifperfile = false;

%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Find list of files

data_dir = fullfile(data_path, 'collated_data/anoved_rsampsl_biprref_cmtspwr');

if isempty(cats)
    fname_cat = 'C';
else
    fname_cat = cats;
end

if isempty(areas)
    fname_area = '*';
else
    fname_area = num2str(areas);
end

fname_out = [fname_cat, '*S' fname_area];

%% Call function
p_thresholder(data_dir, fname_out, q, ifperfile);

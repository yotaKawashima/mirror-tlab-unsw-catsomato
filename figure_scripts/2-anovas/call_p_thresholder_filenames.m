% Calls the p_thresholder_filenames function

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set your selection here
% Set q here
q = 0.05;

% Run per file? false = run for all files together.
ifperfile = false;

%% Find list of files
% find a list of cats (load anova data)
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));
dtypes = {'cmtspwr', 'cmtspwr_snrsurr', 'cmtspwr_evkdpwr'};
for i_dtype = 1:length(dtypes)
    dtype = ['epoched_rsampsl_biprref_evkresp_', dtypes{i_dtype}, ...
             '_adatain_adatout'];
         
    fnames = cell(numel(cat_names)*2, 1);

    for c = 1:numel(cat_names)
        cDir = fullfile(data_path, 'included_datasets', cat_names{c}, dtype);
        files = dir(fullfile(cDir, '*.mat'));
        for k = 1:2
            fnames{(c-1)*2 + k} = fullfile(cDir, files(k).name);
        end
    end % for c = 1:numel(catnames)

    %% Call function
    p_thresholder_filenames(fnames, q, ifperfile);
end % for i_dtype = 1:length(dtypes)

%% Move files into collated_data dir 
%  This part is Added by Yota. 

out_dir = fullfile(data_path, 'collated_data', 'anova_props_test');

% find out if the folder exists.
% mkdir if not exist.
if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end %if ~exist

% move the files
movefile(fullfile('*.mat'), out_dir);
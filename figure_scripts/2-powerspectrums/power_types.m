run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

data_dir = fullfile(data_path, 'included_datasets');
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr';

filename = dir(fullfile('./', ['*' data_type '.mat']));

data = powerplot_data(data_dir, data_type, filename);
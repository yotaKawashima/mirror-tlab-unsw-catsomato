run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

data_dir = fullfile(data_path, 'included_datasets');
premade_data = fullfile(data_path, 'collated_data', 'powerplots_data');

imgformat = '-depsc';

cat_names = dirsinside(data_dir);

%% power
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr';
fname = dir(fullfile(premade_data, ['*' data_type '.mat']));
if ~isempty(fname)
    collectallcats(premade_data, data_type, 'plotmean', 1, 'fromfile', fname.name)
else
    collectallcats(data_dir, data_type, 'plotmean', 1, 'domain', [0, 250])
end
print(1, imgformat, 'Figure2A')

%% SNR
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';
fname = dir(fullfile(premade_data, ['*' data_type '.mat']));
if ~isempty(fname)
    collectallcats(premade_data, data_type, 'plotmean', 2, 'plotttest', 4, 'fromfile', fname.name)
else
    % check if snrsurr files exist
    for c = 1:numel(cat_names)
        loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));

        % check if the data files exist! if not, make them.
        if isempty(loadnames)
            fprint('Creating SNR files\n')
            call_snrtosurrounds(data_path, data_type(1:end-8), cat_names{c})
        end
    end
    
    collectallcats(data_dir, data_type, 'plotmean', 2, 'plotttest', 4)
end
print(2, imgformat, 'Figure2B')
print(4, imgformat, 'Figure2D')

%% EP
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr';
fname = dir(fullfile(premade_data, ['*' data_type '.mat']));
if ~isempty(fname)
    collectallcats(premade_data, data_type, 'plotmean', 3, 'plotttest', 5, 'fromfile', fname.name)
else
    
    % check if files exist
    for c = 1:numel(cat_names)
        loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));

        % check if the data files exist! if not, make them.
        if isempty(loadnames)
            fprint('Creating evoked power files\n')
            call_evokedpower(data_path, data_type(1:end-8), cat_names{c})
        end
    end
    
    collectallcats(data_dir, data_type, 'plotmean', 3, 'plotttest', 5)
end
print(3, imgformat, 'Figure2C')
print(5, imgformat, 'Figure2E')
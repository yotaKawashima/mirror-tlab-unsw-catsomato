%% get overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% set specific variables
% find filenames
data_dir = fullfile(data_path, 'included_datasets'); 
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';
hgp_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp_adatain_adatout_fstonly';

%% Load all of the data 
% Find name for a single cat
cat_names = dirsinside(data_dir);

% preallocate
allcat_fstats = cell(2,1); % holds all of the f-statistics

% frequencies of interest
foi = [23:23:250 200:-23:0];

for c = 1:numel(cat_names)
    % find names
    dir_sec = fullfile(data_dir, cat_names{c}, data_type, '/');
    loadname = dir(fullfile(dir_sec, '*.mat'));
    dir_hgp = fullfile(data_dir, cat_names{c}, hgp_type, '/');
    hgpname = dir(fullfile(dir_hgp, '*.mat'));
    
    for a = 1:2
        % load frequency data
        load(fullfile(dir_sec, loadname(a).name))

        % get fois if needed
        if c == 1 && a == 1
            [~, ind] = find_closest(metavars.freq{1}, foi);
        end
        
        % main fstats
        f_tmp = fstats(:, ind, :);
        
        % load HGP data
        load(fullfile(dir_hgp, loadname(a).name))
        
        % stick it on the end
        f_tmp = [f_tmp, fstats];
        
        % deal with any NaNs
        f_nan = mean(mean(f_tmp, 2), 3); % any NaNs will propogate
        nonnans = find(~isnan(f_nan));
        f_tmp = f_tmp(nonnans, :, :);
        
        % truss
        allcat_fstats{a} = [allcat_fstats{a}; f_tmp];
    end
    
end

%% Plot and save
for a = 1:2
    % plot corr
    figure((a-1)*4+1)
    tmp = horzcat(allcat_fstats{a}(:, :, 1),allcat_fstats{a}(:, :, 2),allcat_fstats{a}(:, :, 3));
    
    tmp1 = corrcoef(tmp);
    
    imagesc(tmp1)
    clim = get(gca, 'CLim');
    ytl = [foi, 0, foi, 0, foi, 0];
    set(gca, 'clim', [-0.06, 1])
    set(gca, 'yticklabel', ytl)
    set(gca, 'ytick', 1:length(ytl))
    set(gca, 'xticklabel', ytl)
    set(gca, 'xtick', 1:length(ytl))
    colorbar

    print(gcf, imgtype, ['correlationmatrix_S' num2str(a)])
    
    % save
    %corrdata = tmp;
    %save(['S' num2str(a) '_corrdata_withNaNs'], 'corrdata')
    
end

%% Save
for a = 1:2
    
end
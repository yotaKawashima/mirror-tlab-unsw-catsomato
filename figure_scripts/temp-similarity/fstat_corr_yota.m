%% get overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% set specific variables
% find filenames
data_dir = fullfile(data_path, 'included_datasets'); 
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';
hgp_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp_adatain_adatout_fstonly';

%imgtype = '-depsc';
imgtype = '-dsvg';

%% Load all of the data 
% Find name for a single cat
cat_names = dirsinside(data_dir);

% preallocate
allcat_fstats = cell(2,1); % holds all of the f-statistics

% find the right frequencies
foi_fund = [23,200];
foi_harm = 46:23:250;
foi_int = sort([177:-23:0 200+23:23:250]);
foi = [foi_fund, foi_harm, foi_int];

n_freq = length(foi) + 1; % # of responses including hgp
n_sessions = numel(cat_names); % # of sessions

% Initialise correlation coeff matrix per session per area.
corr_matrix = zeros(n_freq, n_freq , n_sessions, 3, 2);

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
            [~, ind_fund] = parfind_closest(metavars.freq{1}, foi_fund);
            [~, ind_harm] = parfind_closest(metavars.freq{1}, foi_harm);
            [~, ind_int] = parfind_closest(metavars.freq{1}, foi_int);
            ind = [ind_fund; ind_harm; ind_int];
        end
        
        % main fstats (from frequency of interest)
        % (channels x frequency x 3)
        % 3rd dim for (2 main effects + interaction)?
        f_foi = fstats(:, ind, :);
        
        % load HGP data (fstats from HGP)
        % (1 x channels x 3) 
        % 1st dim for HGP : HGP is defined as mean within a frequency
        % range
        % 3rd dim for (2 main effects + interaction)?
        load(fullfile(dir_hgp, hgpname(a).name))
        % convert it into (channels x HGP x 3).
        f_hgp = permute(fstats, [2, 1, 3]);
        
        % combine two fstats
        f_all_tmp = [f_foi, f_hgp];
        
        % deal with any NaNs
        % take mean across frequency+HGP across 3 to remove channels which
        % generate nan value.
        f_nan = mean(mean(f_all_tmp, 2), 3); % any NaNs will propogate
        nonnans = find(~isnan(f_nan)); 
        f_all = f_all_tmp(nonnans, :, :); % remove nan channels 
        
        % convert it into channels x ((foi + hgp)x3)
        %corrdata = horzcat(f_all(:, :, 1), f_all(:, :, 2), f_all(:, :, 3));
        
        % compute correlation (why take log of f-statistic?)
        %[tmpl, ~] = corrcoef(log(corrdata));
        %corr_matrix(:,:,c,a) = tmpl;

        for i = 1:3
            corrdata = f_all(:,:,i);
            [tmpl, ~] = corrcoef(log(corrdata));
            %tmpl(1:1+size(tmpl,1):end) = nan;
            corr_matrix(:,:,c,i,a) = tmpl;
        end % i
    end % a
end % c

%% Plot mean correlation coeff matrix per area
for a = 1:2
    for i = 1:3
        % take mean across sessions
        %mean_corrmatrix = mean(corr_matrix(:,:,:,a) , 3);
        mean_corrmatrix = mean(corr_matrix(:,:,:,i,a), 3);
        figure(a)

        imagesc(mean_corrmatrix)
        clim = get(gca, 'CLim');
        ytl = [foi, 0]; % 0 is HGP
        set(gca, 'clim', [0, 1])
        set(gca, 'yticklabel', ytl)
        set(gca, 'ytick', 1:length(ytl))
        set(gca, 'xticklabel', ytl)
        set(gca, 'xtick', 1:length(ytl))
        colorbar

        print(gcf, imgtype, ['mean_correlationmatrix_S' num2str(a), '_fid_', num2str(i)])

        % save
        %save(['S' num2str(a) '_corrdata_withNaNs'], 'corrdata')
    end
end

%% Plot mean correlation coeff matrix collapsing harmonic and intermodulation per area
for a = 1:2
    for i = 1:3
        ytl = ['f1'; 'f2'; 'f1harmonic'; 'intermodulation'; 'HGB']; 

        mean_collapse_corrmatrix = zeros(length(ytl), length(ytl));
        
        % take mean across sessions
        mean_corrmatrix = mean(corr_matrix(:,:,:,i,a) , 3);
        
        mean_collapse_corrmatrix(1,1) = mean_corrmatrix;
        mean_collapse_corrmatrix(1,2) = mean_corrmatrix(1,2);
        mean_collapse_corrmatrix(1,3) = mean_corrma
            
        figure(a)

        imagesc(mean_corrmatrix)
        clim = get(gca, 'CLim');
        set(gca, 'clim', [0, 1])
        set(gca, 'yticklabel', ytl)
        set(gca, 'ytick', 1:length(ytl))
        set(gca, 'xticklabel', ytl)
        set(gca, 'xtick', 1:length(ytl))
        colorbar

        print(gcf, imgtype, ['mean_collapse_correlationmatrix_S' num2str(a), '_fid_', num2str(i)])

        % save
        %save(['S' num2str(a) '_corrdata_withNaNs'], 'corrdata')
    end
end
    

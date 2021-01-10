%% Get top 10% channels for model fitting.
% Fing top 10% of channels that give larger value of 
% sum(logSNR at harmonics and IMs).

%% Set overall variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));

% select datatype
datatype_logsnr = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

%img_fmt = '-dpng';
img_fmt = '-depsc';

%% Harmonics and Intermodulations
f1 = 23;
f2 = 200;
f1_all = 23:23:250;
f1_harms = 46:23:250;
ims = sort([177:-23:0 200+23:23:250]);
harm_im = sort([f1_harms , ims]);     
foi_all = sort([f1_all, f2, ims]);            

%% Get logSNR from all channels (= bipolar channels x sessions) for each area
% Get only the max vibration condition from each session.
% Initialisation
data_struct = struct();

for area_id = 1:2
    area = ['S', num2str(area_id)]; % area
    
    % Initialisation 
    mean_logsnr_0_250_matrix = [];
    sd_logsnr_0_250_matrix = [];
    mean_logsnr_at_harms_ims_matrix = [];
    session_matrix = [];
    bipolar_chs_matrix =[];
    
    for session_id = 1:numel(cat_names)
        % find file names
        session_path = fullfile(data_dir, cat_names{session_id}, datatype_logsnr);
        loadname = dir(fullfile(session_path, ['*', area, '*.mat']));
        % load only the max vibration condition
        load(fullfile(session_path, loadname(end).name));
        
        % Get logSNR data
        % check if unipolar channels are included. If so, remove the
        % unipolar channels.
        if strcmp(data.label{1}(1:3), 'raw')
            % the first bipolar channel
            bipolarch_beginning = 1 + prod(data.custom.spatialconfig);
            % # of bipolar channels
            nChan = data.custom.nsignals - bipolarch_beginning + 1; 
        else
            bipolarch_beginning = 1;
            nChan = data.custom.nsignals;
        end
        
        % Get only bipolar channels data.
        % data dim = (bipolar channels x freqs x trials)
        logsnrs_ = data.trial(bipolarch_beginning:end, :, :);
        
        % Get logsnrs from 0Hz to 250Hz
        [~, fmax_ind] = find_closest(data.freq{1}, 250);
        logsnrs_0_250 = logsnrs_(:, 1:fmax_ind, :);
        
        % Mean across trials
        % data dim = bipolar channels x freqs
        logsnrs_0_250_mean_across_trials = squeeze(mean(logsnrs_0_250, 3));
        logsnrs_0_250_sd_across_trials = squeeze(std(logsnrs_0_250, 0, 3));
        
        % Get freq ids corresponding to harmonics and IMs
        [~, harm_im_inds] = find_closest(data.freq{1}, harm_im);
        
        % Get logSNR at harmonics and IMs
        % data dim = bipolar channels x (harmonics and ims)
        logsnrs_at_harm_im_mean_across_trials = ...
                logsnrs_0_250_mean_across_trials(:, harm_im_inds, :);
        
        % Store data across sessions
        % data dim = (bipolar channels x sessions) x (harmonics and ims)
        mean_logsnr_at_harms_ims_matrix = [mean_logsnr_at_harms_ims_matrix; ...
                                      logsnrs_at_harm_im_mean_across_trials];
        % data dim = (bipolar channels x sessions) x freqs
        mean_logsnr_0_250_matrix = [mean_logsnr_0_250_matrix;...
                                    logsnrs_0_250_mean_across_trials]; 
        sd_logsnr_0_250_matrix = [sd_logsnr_0_250_matrix;...
                                  logsnrs_0_250_sd_across_trials];
        session_matrix = [session_matrix;...
                          repmat(cat_names{session_id}, nChan, 1)];
        bipolar_chs_matrix = [bipolar_chs_matrix; (1:nChan)'];
    end % for session_id = 1:numel(cat_names)
    % Here, take sum of logsnrs across harms and ims. Based on the sum,
    % take top ??% channels.
    sum_logsnrs_matrix = sum(mean_logsnr_at_harms_ims_matrix, 2);
    
    % Store data with structure for each area
    data_struct(area_id).area = area;
    data_struct(area_id).freqs = data.freq{1}(1:fmax_ind);
    data_struct(area_id).session = session_matrix;
    data_struct(area_id).bipolar_ch = bipolar_chs_matrix;
    data_struct(area_id).sum_logsnrs_matrix = sum_logsnrs_matrix;
    data_struct(area_id).mean_logsnr_at_harms_ims_matrix = mean_logsnr_at_harms_ims_matrix;
    data_struct(area_id).mean_logsnr_0_250_matrix = mean_logsnr_0_250_matrix;
    data_struct(area_id).sd_logsnr_0_250_matrix = sd_logsnr_0_250_matrix;
    
end % for area_id = 1:2


%% Get thresholds for each area
% top 5% 
top10_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 90, 1);
top10_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 90, 1);

% Get top 5% channels for each area 
top10_channels_S1 = find(data_struct(1).sum_logsnrs_matrix > top10_threshold_S1);
top10_channels_S2 = find(data_struct(2).sum_logsnrs_matrix > top10_threshold_S2);

% Get log SNR of the 5% (mean across trials)
top10_mean_logsnrs_S1 = data_struct(1).mean_logsnr_0_250_matrix(top10_channels_S1, :);
top10_mean_logsnrs_S2 = data_struct(2).mean_logsnr_0_250_matrix(top10_channels_S2, :);

% Store data for each area
data_top10_struct = struct();
data_top10_struct(1).top_channels = top10_channels_S1;
data_top10_struct(2).top_channels = top10_channels_S2;
data_top10_struct(1).session = data_struct(1).session(top10_channels_S1, :);
data_top10_struct(2).session = data_struct(2).session(top10_channels_S2, :);
data_top10_struct(1).bipolar_ch = data_struct(1).bipolar_ch(top10_channels_S1, 1);
data_top10_struct(2).bipolar_ch = data_struct(2).bipolar_ch(top10_channels_S2, 1);
data_top10_struct(1).top_threshold = top10_threshold_S1;
data_top10_struct(2).top_threshold = top10_threshold_S2;
data_top10_struct(1).mean_logsnrs = top10_mean_logsnrs_S1;
data_top10_struct(2).mean_logsnrs = top10_mean_logsnrs_S2;

% Save data
save('Harmonics_IMs_top10.mat', 'data_top10_struct');

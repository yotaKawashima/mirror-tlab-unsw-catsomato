%% Check harmonics and IMs channels
% Find top 10% of channels that give the larger value of sum(logSNR at
% harmonics and IMs).


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
f1_all = 23:23:250;
f1_harms = 46:23:250;
ims = sort([177:-23:0 200+23:23:250]);
harm_im = sort([f1_harms , ims]);     

f2 = 200;
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

    end % for session_id = 1:numel(cat_names)
    % Here, take sum of logsnrs across harms and ims. Based on the sum,
    % take top ??% channels.
    sum_logsnrs_matrix = sum(mean_logsnr_at_harms_ims_matrix, 2);
    
    % Store data with structure for each area
    data_struct(area_id).area = area;
    data_struct(area_id).sum_logsnrs_matrix = sum_logsnrs_matrix;
    data_struct(area_id).mean_logsnr_at_harms_ims_matrix = mean_logsnr_at_harms_ims_matrix;
    data_struct(area_id).mean_logsnr_0_250_matrix = mean_logsnr_0_250_matrix;
    data_struct(area_id).sd_logsnr_0_250_matrix = sd_logsnr_0_250_matrix;
    
end % for area_id = 1:2

%% Plot thresholds for each area
% top 5% 
top5_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 95, 1);
top5_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 95, 1);

% top 10%
top10_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 90, 1);
top10_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 90, 1);

% top 25%
top25_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 75, 1);
top25_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 75, 1);

% top 50%
top50_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 50, 1);
top50_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 50, 1);

S1_thresholds =  [top5_threshold_S1, top10_threshold_S1,...
                  top25_threshold_S1, top50_threshold_S1];
              
S2_thresholds = [top5_threshold_S2, top10_threshold_S2, ...
                 top25_threshold_S2, top50_threshold_S2];
figure();clf              
hold on 
scatter([0.25, 0.75], [S1_thresholds(1), S2_thresholds(1)], 100, 'filled');
scatter([0.25, 0.75], [S1_thresholds(2), S2_thresholds(2)], 100, 'filled');
scatter([0.25, 0.75], [S1_thresholds(3), S2_thresholds(3)], 100, 'filled');
scatter([0.25, 0.75], [S1_thresholds(4), S2_thresholds(4)], 100, 'filled');
xlim([0, 1]);
xticks([0.25, 0.75]);
xticklabels({'S1', 'S2'});
ylabel({'sum of logSNR';'across harmonics and IMs'});
set(gca, 'fontsize', 16);
legend({'top 5% threshold', 'top 10% threshold', ...
        'top 25% threshold', 'top 50% threshold'}, 'location', 'northeastoutside');
hold off


%% Plot top 5% channels for each area
% Get top 5% channels for each area 
top5_channels_S1 = find(data_struct(1).sum_logsnrs_matrix > top5_threshold_S1);
top5_channels_S2 = find(data_struct(2).sum_logsnrs_matrix > top5_threshold_S2);

% Frequency data
freq_0_250 = data.freq{1}(1:fmax_ind);

% color for shading
col = [192,192,192]/256;

% Plot S1 logSNR first
top5_mean_logsnrs_S1 = data_struct(1).mean_logsnr_0_250_matrix(top5_channels_S1, :);
top5_sd_logsnrs_S1 = data_struct(1).sd_logsnr_0_250_matrix(top5_channels_S1, :);

increment = 1;
increment_fig = 0;
for ch_id = 1:length(top5_channels_S1)
    % Create new figure 
    if or(ch_id == 1, increment > 20)
        figure();clf
        increment_fig = increment_fig + 1;
        title_ = ['fig ID ', num2str(increment_fig), ...
                  ', S1 top 5%: ', num2str(length(top5_channels_S1)), '/', ...
                  num2str(size(data_struct(1).mean_logsnr_0_250_matrix, 1)),...
                  'chs'];
        sgtitle(title_);
        increment = 1; % update increment for subplot
    end
    
    subplot(5, 4, increment)
    hold on
    up_ = top5_mean_logsnrs_S1(ch_id, :) + top5_sd_logsnrs_S1(ch_id, :);
    low_ = top5_mean_logsnrs_S1(ch_id, :) - top5_sd_logsnrs_S1(ch_id, :);    
    h=fill([freq_0_250', fliplr(freq_0_250')], [up_, fliplr(low_)], 'r');
    set(h, 'FaceColor', col)
    set(h, 'EdgeColor', col)
    plot(freq_0_250, top5_mean_logsnrs_S1(ch_id, :), 'k');
    xlim([0, 250]);
    ylim([-5, 20]);   
    
    % plot vertical lines
    v_ylims = [-5, -3];
    v1 = numel(f1_all);
    v2 = numel(ims);
    v = zeros(1, numel(foi_all));
    v(1:v1) = plot([f1_all', f1_all'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
    v(v1+1) = plot([f2, f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
    v(v1+2:v1+v2+1) = plot([sort(ims)', sort(ims)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
    % format vertical lines
    set(v, 'LineWidth', 2, 'LineStyle', '--')    
    hold off
    
    increment = increment + 1;
            
end % ch_id =1:length(top_channels_S1)



% Plot S2 logSNR second
top5_mean_logsnrs_S2 = data_struct(2).mean_logsnr_0_250_matrix(top5_channels_S2, :);
top5_sd_logsnrs_S2 = data_struct(2).sd_logsnr_0_250_matrix(top5_channels_S2, :);

increment = 1;
increment_fig = 0;
for ch_id = 1:length(top5_channels_S2)
    % Create new figure 
    if or(ch_id == 1, increment > 20)
        figure();clf
        increment_fig = increment_fig + 1;
        title_ = ['fig ID ', num2str(increment_fig), ...
                  ', S2 top 5%: ', num2str(length(top5_channels_S2)), '/', ...
                  num2str(size(data_struct(2).mean_logsnr_0_250_matrix, 1)),...
                  'chs'];
        sgtitle(title_);
        increment = 1; % update increment for subplot
    end
    
    subplot(5, 4, increment)
    hold on
    up_ = top5_mean_logsnrs_S2(ch_id, :) + top5_sd_logsnrs_S2(ch_id, :);
    low_ = top5_mean_logsnrs_S2(ch_id, :) - top5_sd_logsnrs_S2(ch_id, :);    
    h=fill([freq_0_250', fliplr(freq_0_250')], [up_, fliplr(low_)], 'r');
    set(h, 'FaceColor', col)
    set(h, 'EdgeColor', col)
    plot(freq_0_250, top5_mean_logsnrs_S2(ch_id, :), 'k');
    xlim([0, 250]);
    ylim([-5, 20]);   
    
    % plot vertical lines
    v_ylims = [-5, -3];
    v1 = numel(f1_all);
    v2 = numel(ims);
    v = zeros(1, numel(foi_all));
    v(1:v1) = plot([f1_all', f1_all'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
    v(v1+1) = plot([f2, f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
    v(v1+2:v1+v2+1) = plot([sort(ims)', sort(ims)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
    % format vertical lines
    set(v, 'LineWidth', 2, 'LineStyle', '--')    
    hold off
    
    increment = increment + 1;                
    
end % ch_id =1:length(top_channels_S2)




      
%% Functions

%{
    % Define the initial simplex as tabel
    coordinate = [0, 0; var1max, 0; 0, var2max];
    error1 = inf;
    error2 = compute_error(model(var1max, 0), data, sampling_rate);
    error3 = compute_error(model(0, var2max), data, sampling_rate); 
    error = [error1; error2; error3];
    simplex = table(coordinate, error);
    
    % Get the best vertix from the last simplex
    simplex = sortrows(simplex, 'error');  
    min_var1 = simplex.coordinate(1, 1); 
    min_var2 = simplex.coordinate(1, 2);
    error = simplex.error(1);
%}
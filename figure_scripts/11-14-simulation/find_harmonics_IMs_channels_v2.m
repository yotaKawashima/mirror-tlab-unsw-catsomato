%% Check harmonics and IMs channels
% Find top 10% of channels that give the larger value of sum(logSNR at
% harmonics and IMs).
% Plot histogram of each responses.

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

%% Get thresholds for each area
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

% Get top 5% channels for each area 
top5_channels_S1 = find(data_struct(1).sum_logsnrs_matrix > top5_threshold_S1);
top5_channels_S2 = find(data_struct(2).sum_logsnrs_matrix > top5_threshold_S2);

% Get log SNR of the 5% (mean across trials)
top5_mean_logsnrs_S1 = data_struct(1).mean_logsnr_0_250_matrix(top5_channels_S1, :);
top5_mean_logsnrs_S2 = data_struct(2).mean_logsnr_0_250_matrix(top5_channels_S2, :);

%% Histogram of ratio (samples --- top 5% channels)
% Get ratio 23 vs 46 
[~, f23_ind] = find_closest(data.freq{1}, 23);
[~, f46_ind] = find_closest(data.freq{1}, 46);
ratio_23_46_S1 = top5_mean_logsnrs_S1(:, f23_ind) ... 
                 ./ top5_mean_logsnrs_S1(:, f46_ind);
ratio_23_46_S2 = top5_mean_logsnrs_S2(:, f23_ind) ... 
                 ./ top5_mean_logsnrs_S2(:, f46_ind);

% Get ratio 23 vs 69 
[~, f69_ind] = find_closest(data.freq{1}, 69);
ratio_23_69_S1 = top5_mean_logsnrs_S1(:, f23_ind) ...
                 ./ top5_mean_logsnrs_S1(:, f69_ind);
ratio_23_69_S2 = top5_mean_logsnrs_S2(:, f23_ind) ...
                 ./ top5_mean_logsnrs_S2(:, f69_ind);

% Get ratio 200 vs 177 
[~, f200_ind] = find_closest(data.freq{1}, 200);
[~, f177_ind] = find_closest(data.freq{1}, 177);
ratio_200_177_S1 = top5_mean_logsnrs_S1(:, f200_ind) ...
                 ./ top5_mean_logsnrs_S1(:, f177_ind);
ratio_200_177_S2 = top5_mean_logsnrs_S2(:, f200_ind) ... 
                 ./ top5_mean_logsnrs_S2(:, f177_ind);
%{
figure(); clf
subplot(2, 3, 1)
histogram(ratio_23_46_S1);
title('S1 response 23Hz / response 46HZ');
subplot(2, 3, 2)
histogram(ratio_23_69_S1);
title('S1 response 23Hz / response 69HZ');
subplot(2, 3, 3)
histogram(ratio_200_177_S1);
title('S1 response 200Hz / response 177HZ');
subplot(2, 3, 4)
histogram(ratio_23_46_S2);
title('S2 response 23Hz / response 46HZ');
subplot(2, 3, 5)
histogram(ratio_23_46_S2);
title('S2 response 23Hz / response 69HZ');
subplot(2, 3, 6)
histogram(ratio_200_177_S2);
title('S2 response 200Hz / response 177HZ');
%}
                 
%{
figure(); clf
subplot(1, 3, 1)
hold on
h1 = histogram(ratio_23_46_S1);
h2 = histogram(ratio_23_46_S2);
h1.Normalization = 'probability';
h1.BinWidth = 1;
h2.Normalization = 'probability';
h2.BinWidth = 1;
hold off
title('logSNR 23Hz / logSNR 46HZ');
subplot(1, 3, 2)
hold on
h1 = histogram(ratio_23_69_S1);
h2 = histogram(ratio_23_69_S2);
h1.Normalization = 'probability';
h1.BinWidth = 1;
h2.Normalization = 'probability';
h2.BinWidth = 1;
hold off
title('log SNR 23Hz / logSNR 69HZ');
subplot(1, 3, 3)
hold on
h1 = histogram(ratio_200_177_S1);
h2 = histogram(ratio_200_177_S2);
h1.Normalization = 'probability';
h1.BinWidth = 1;
h2.Normalization = 'probability';
h2.BinWidth = 1;
hold off
title('logSNR 200Hz / logSNR 177HZ');
legend({'S1', 'S2'});
%}
                 
figure(); clf
subplot(1, 3, 1)
hold on
h1 = histogram(ratio_23_46_S1);
h2 = histogram(ratio_23_46_S2);
h1.BinWidth = 1;
h2.BinWidth = 1;
hold off
title('logSNR 23Hz / logSNR 46HZ');
subplot(1, 3, 2)
hold on
h1 = histogram(ratio_23_69_S1);
h2 = histogram(ratio_23_69_S2);
h1.BinWidth = 1;
h2.BinWidth = 1;
hold off
title('log SNR 23Hz / logSNR 69HZ');
subplot(1, 3, 3)
hold on
h1 = histogram(ratio_200_177_S1);
h2 = histogram(ratio_200_177_S2);
h1.BinWidth = 1;
h2.BinWidth = 1;
title('logSNR 200Hz / logSNR 177HZ');
legend({'S1', 'S2'});

%% Plot scatter
figure(); clf
hold on
scatter(ratio_23_46_S1, ratio_23_69_S1, 'filled');
scatter(ratio_23_46_S2, ratio_23_69_S2, 'filled');
hold off 
xlabel('logSNR 23Hz / logSNR 46Hz');
ylabel('logSNR 23Hz / logSNR 69Hz');
legend({'S1', 'S2'});
set(gca, 'fontsize', 16);

figure(); clf
hold on
scatter(ratio_23_46_S1, ratio_200_177_S1, 'filled');
scatter(ratio_23_46_S2, ratio_200_177_S2, 'filled');
hold off
xlabel('logSNR 23Hz / logSNR 46Hz');
ylabel('logSNR 200Hz / logSNR 177Hz');
legend({'S1', 'S2'});
set(gca, 'fontsize', 16);

figure(); clf
hold on
scatter(ratio_23_69_S1, ratio_200_177_S1, 'filled');
scatter(ratio_23_69_S2, ratio_200_177_S2, 'filled');
hold off
xlabel('logSNR 23Hz / logSNR 69Hz');
ylabel('logSNR 200Hz / logSNR 177Hz');
legend({'S1', 'S2'});
set(gca, 'fontsize', 16);


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
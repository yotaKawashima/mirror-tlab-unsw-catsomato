% Plot VELogP for a wider range. 
% Plots the evoked power in each condition across 50-150 Hz
% The condition (F1 amplitude ,F2 amplitude) = (0, 0) is included.
%% Set variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

cat_name = 'C20110808_R03';
bipolar_channel_id = 102;
x_width = 10;
y_width = 6;
ftype = 'Arial';
fsize = 10;
img_fmt = '-dpdf';

% Data directory
data_dir = fullfile(data_path, 'included_datasets');

%% Get top 10 % channels (consider only HGP range)
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp';

% find names
cat_names = dirsinside(data_dir);
loadname = dir(fullfile(data_dir, cat_names{1}, datatype, '*_S2_*.mat'));
for k = 2:numel(cat_names)
    loadname(k) = ...
        dir(fullfile(data_dir, cat_names{k}, datatype, '*_S2_*.mat'));
end

% change to store hgp
mean_allhgp = [];
sessions = [];
bpchannel_matrix = [];
% load HGP data from sessions and store it.
for k = 1:numel(loadname)
    load(fullfile(data_dir, cat_names{k}, datatype, loadname(k).name))

    % work out which channels to keep (the bipolar ones)
    if data.custom.nsignals == 280 || data.custom.nsignals == 176
        lastunipol = prod(data.custom.spatialconfig);
        chan = lastunipol+1:data.custom.nsignals;
    elseif data.custom.nsignals == 180 || data.custom.nsignals == 112
        chan = 1:data.custom.nsignals;
    else
        error('Only unipolar?');
    end
    % Get the maximal stimulus condition
    max_data = data.trial{end};
    
    % Compute HGP i.e. mean across frequencies 
    % (from 50Hz t0 150Hz excluding foi).
    % chan x freqs x trials
    max_hgp_data = squeeze(mean(max_data, 2));
    
    % Extract only bipolar channel data
    % now bp chan x trials
    max_hgp_data_bp = max_hgp_data(chan, :);
    
    % Mean/std/stderr across trials per bp chan
    mean_max_hgp_data_bp = mean(max_hgp_data_bp, 2);
    
    mean_allhgp = [mean_allhgp; mean_max_hgp_data_bp]; % bp chans 
    sessions = [sessions ; ones(length(chan),1)*k];
    bpchannel_matrix = [bpchannel_matrix, 1:length(chan)];
end % for k = 1:numel(loadname)

% Get top 10 channels (remove nan chans)
mean_hgp_without_nan = mean_allhgp(~isnan(mean_allhgp));
prct90 = prctile(mean_hgp_without_nan, 90, 1);
top_hgp_ch = mean_allhgp >= prct90;

%% Collect logpower data ans compute VELogP
% load logpower
datatype_logpow = 'epoched_rsampsl_biprref_evkresp_cmtspwr';

all_mean_velogp = [];
all_stderr_velogp = [];

for cat_id = 1:numel(cat_names)
    cat_name = cat_names{cat_id}; 
    % find names (list files correspoing to each stim condition)
    loadnames = dir(fullfile(data_dir, cat_name, datatype_logpow, '*.mat'));

    % Get the min and the max stimulus vibration condition in S2. 
    % (the first and the last file)
    min_condid = numel(loadnames)/2 + 1;
    max_condid = numel(loadnames);
    min_stim_data = ...
        load(fullfile(data_dir, cat_name, datatype_logpow, loadnames(min_condid).name));
    max_stim_data = ...
        load(fullfile(data_dir, cat_name, datatype_logpow, loadnames(max_condid).name));

    % check if unipolar data is included or not.
    if strcmp(min_stim_data.data.label{1}(1:3), 'raw')
        % the first bipolar channel
        bipolarchid = 1 + prod(min_stim_data.data.custom.spatialconfig);
        % # of bipolar channels
        nChan = min_stim_data.data.custom.nsignals - bipolarchid + 1; 
    else
        nChan = data.custom.nsignals;
        bipolarchid = 1;
    end
    % Extract responses of bipolar channels
    % chans x freqs x trials
    min_stim_data_bp = min_stim_data.data.trial(bipolarchid:end, :, :);
    max_stim_data_bp = max_stim_data.data.trial(bipolarchid:end, :, :);

    % compute baseline (mean accross trials in no stim cond)
    baseline = mean(min_stim_data_bp, 3);

    % Compute VELogP for max stim cond
    velogp_max_stim_bp = max_stim_data_bp - baseline;
    
    % Mean/std/stderr across trials
    % channels x freqs x trials
    mean_velogp_max_stim_bp = squeeze(mean(velogp_max_stim_bp, 3));
    std_velogp_max_stim_bp = squeeze(std(velogp_max_stim_bp, 0, 3));
    stderr_velogp_max_stim_bp = ...
        std_velogp_max_stim_bp / sqrt(size(velogp_max_stim_bp, 3));
    
    % Store data
    % channels x freqs
    all_mean_velogp = ...
        [all_mean_velogp; ...
         mean_velogp_max_stim_bp];
    all_stderr_velogp = ...
        [all_stderr_velogp;...
         stderr_velogp_max_stim_bp];
end % cat_id = 1:numel(cat_names)

% Get Nyqyuist freq.
Ny_freq = min_stim_data.data.fsample/2;

% frequencies
freqs = squeeze(min_stim_data.data.freq{1});

% Get data of top 10% channels
top10_mean_velogp = all_mean_velogp(top_hgp_ch, :);
top10_stderr_velogp = all_stderr_velogp(top_hgp_ch, :);

% Mean across top 10 %channels
mean_mean_velogp = mean(top10_mean_velogp);
mean_stderr_velogp = mean(top10_stderr_velogp);

%% Plot figures
% Remove 50Hz noise
multiple_50 = [1:1:99] * 50 ;
[~, multiple_50_ids] = ...
    find_closest(freqs, multiple_50);

% remove harmonics
multiple_23 = [23:23:Ny_freq];
[~, multiple_23_ids] = ...
    find_closest(freqs, multiple_23);

% remove some intermodulats
intermoduls = sort([200:-23:0 200:23:4000]);
[~, intermoduls_ids] = ...
    find_closest(freqs, intermoduls);

% Replace XHz+-~3Hz multiple with nan 
% freq step = 0.15Hz 
for freq_id = -10:10
    mean_mean_velogp(multiple_50_ids+freq_id) = nan;  
    mean_mean_velogp(multiple_23_ids+freq_id) = nan;
    mean_mean_velogp(intermoduls_ids+freq_id) = nan;
    
    mean_stderr_velogp(multiple_50_ids+freq_id) = nan;  
    mean_stderr_velogp(multiple_23_ids+freq_id) = nan;
    mean_stderr_velogp(intermoduls_ids+freq_id) = nan;
    
end % freq_id = -20:20

figure();clf
hold on
upper = squeeze(mean_mean_velogp + mean_stderr_velogp);
lower = squeeze(mean_mean_velogp - mean_stderr_velogp);
upper = fillmissing(upper,'linear');% 'movmedian', 30);
lower = fillmissing(lower, 'linear');%     'movmedian', 30);

hold on
% Show std with gray shading
h = fill([freqs', fliplr(freqs')], ...
    [upper, fliplr(lower)], 'r');
set(h, 'FaceColor', [192,192,192]/256);
set(h, 'EdgeColor', [192,192,192]/256);
% plot mean
plot(freqs, mean_mean_velogp, 'k');
% plot y=0 line
yline(0, '-k', 'linewidth', 0.5);
xline(50, '-r', 'linewidth', 1);
xline(150, '-r', 'linewidth', 1);
hold off
xlim_upper = Ny_freq; 
xlim([0 xlim_upper]);
xlabel('frequency [Hz]');
ylabel('VELogP [dB]')

set(gca, 'FontSize', fsize);
set(gca, 'FontName', ftype);
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc" || img_fmt == "-dpdf"  || img_fmt == "-dsvg" 
    print(f, img_fmt, 'FigRev3_VELogP_top10');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_VELogP_top10', '-r400');
end

% log narrow
xlim([0 1000]);
xticks([0:100:1000]);
if img_fmt == "-depsc" || img_fmt == "-dpdf" || img_fmt == "-dsvg"
    print(f, img_fmt, 'FigRev3_VELogP_top10_narrow');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_VELogP_top10_narrow', '-r400');
end

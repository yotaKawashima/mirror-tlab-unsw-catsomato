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
y_width = 8;
ftype = 'Arial';
fsize = 10;
img_fmt = '-dpng';

% Data directory
data_dir = fullfile(data_path, 'included_datasets');
%% First extract data and compute VELogP
% load logpower
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr';
loadnames = dir(fullfile(data_dir, cat_name, datatype, '*.mat'));
% Get the min and the max stimulus vibration condition in S2. 
% (the first and the last file)
min_condid = numel(loadnames)/2 + 1;
max_condid = numel(loadnames);
min_stim_data = ...
    load(fullfile(data_dir, cat_name, datatype, loadnames(min_condid).name));
max_stim_data = ...
    load(fullfile(data_dir, cat_name, datatype, loadnames(max_condid).name));

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
min_stim_data_bp = min_stim_data.data.trial(bipolarchid:end, :, :);
max_stim_data_bp = max_stim_data.data.trial(bipolarchid:end, :, :);

% Target channel 
min_stim_data_bp_onech = ...
    squeeze(min_stim_data_bp(bipolar_channel_id, :, :));
max_stim_data_bp_onech = ...
    squeeze(max_stim_data_bp(bipolar_channel_id, :, :));

% compute baseline (mean accross trials in no stim cond)
baseline = mean(min_stim_data_bp_onech, 2);

% Compute VELogP for max stim cond
velogp_max_stim_onech = max_stim_data_bp_onech - baseline;

% Get Nyqyuist freq.
Ny_freq = min_stim_data.data.fsample/2;

% frequencies
freqs = min_stim_data.data.freq{1};

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
intermoduls = sort([577:-23:0 200+23:23:500]);
[~, intermoduls_ids] = ...
    find_closest(freqs, intermoduls);

% Replace 50Hz+-~3Hz multiple with nan 
% freq step = 0.15Hz 
for freq_id = -10:10
    velogp_max_stim_onech(multiple_50_ids+freq_id, :) = nan;  
    velogp_max_stim_onech(multiple_23_ids+freq_id, :) = nan;
    velogp_max_stim_onech(intermoduls_ids+freq_id, :) = nan;
    
end % freq_id = -20:20

% mean/std across trials
mean_velogp = nanmean(velogp_max_stim_onech, 2);
std_velogp = nanstd(velogp_max_stim_onech, 0, 2);
%stderr_velogp = std_velogp / sqrt(size(velogp_max_stim_onech, 2));

figure();clf
hold on
upper = mean_velogp + std_velogp;
lower = mean_velogp - std_velogp;
upper = fillmissing(upper,'linear');% 'movmedian', 30);
lower = fillmissing(lower, 'linear');%     'movmedian', 30);

hold on
% Show std with gray shading
h = fill([freqs', fliplr(freqs')], ...
    [upper', fliplr(lower')], 'r');
set(h, 'FaceColor', [192,192,192]/256);
set(h, 'EdgeColor', [192,192,192]/256);
% plot mean
plot(freqs, mean_velogp, 'k');
% plot y=0 line
yline(0, '-w', 'linewidth', 2);
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
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(f, img_fmt, 'FigRev3_VELogP');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_VELogP', '-r300');
end

% log narrow
xlim([0 1000]);
xticks([0:100:1000]);
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(f, img_fmt, 'FigRev3_VELogP_narrow');moduls
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_VELogP_narrow', '-r300');
end

% log Scale
%xlim([0 xlim_upper]);
%set(gca, 'XScale', 'log')
%if img_fmt == "-depsc" || img_fmt == "-dpdf"   
%    print(f, img_fmt, 'FigRev3_VELogP_log');
%elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
%    print(f, img_fmt, 'FigRev3_VELogP_log', '-r300');
%end
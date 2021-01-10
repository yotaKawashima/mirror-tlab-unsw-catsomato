% Plot power spectrums for pre-stimulus and post-stimulus
% See differences in HGP between them.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

cat_name = 'C20110808_R03';
bipolar_channel_id = 102;
img_fmt = '-dpng';

% Data directory
data_dir = fullfile(data_path, 'included_datasets');
%% First extract pre-stimulus data and compute power spectrum
% load bipolar data
datatype_biprref = 'epoched_rsampsl_biprref';
loadnames = dir(fullfile(data_dir, cat_name, datatype_biprref, '*.mat'));
% Get the max stimulus vibration condition in S2. (the last file)
condid = numel(loadnames); 
load(fullfile(data_dir, cat_name, datatype_biprref, loadnames(condid).name));

% Extract the data from -0.7s to 0s.
prestim_ids = find(data.time{1} <= 0);
sampling_rate = data.fsample;
% check if unipolar data is included or not.
if strcmp(data.label{1}(1:3), 'raw')
    % the first bipolar channel
    bipolarchid = 1 + prod(data.custom.spatialconfig);
    % # of bipolar channels
    nChan = data.custom.nsignals - bipolarchid + 1; 
else
    nChan = data.custom.nsignals;
    bipolarchid = 1;
end
data_biprref = data.trial(bipolarchid:end, :, :);
% time x trials
data_biprref_prestim = ...
    squeeze(data_biprref(bipolar_channel_id, 1:prestim_ids(end), :));

% Compute power spectrum for the pre-stimulus data.
[pow_spect_prestim, freqs_prestim] = ...
    power_spect(data_biprref_prestim, data.fsample);

% Get Nyqyuist freq.
Ny_freq_prestim = data.fsample/2;

%% Get post-stimulus (from 0.5s to 2.5s) power spectrum.
% The power spectrum has already been stored. 
datatype_cmtsgrm = 'epoched_rsampsl_biprref_evkresp_cmtspwr';
    
% Find files
loadnames = dir(fullfile(data_dir, cat_name, datatype_cmtsgrm, '*.mat'));
% Get the max stimulus vibration condition in S2. (the last file)
condid = numel(loadnames); 
load(fullfile(data_dir, cat_name, datatype_cmtsgrm, loadnames(condid).name));

% Get the power spectrum we want.
% check if unipolar data is included or not.
if strcmp(data.label{1}(1:3), 'raw')
    % the first bipolar channel
    bipolarchid = 1 + prod(data.custom.spatialconfig);
    % # of bipolar channels
    nChan = data.custom.nsignals - bipolarchid + 1; 
else
    nChan = data.custom.nsignals;
    bipolarchid = 1;
end
data_bipolar = data.trial(bipolarchid:end, :, :);
% time x trials
pow_spect_pststim= ...
    squeeze(data_bipolar(bipolar_channel_id, :, :));
freqs_pststim = data.freq{1};
% Get Nyquist freq
Ny_freq_pststim = data.fsample/2;

%% Plot figures
% Remove 50Hz noise
multiple_50 = [1:1:99] * 50 ;
[~, multiple_50_ids_prestim] = ...
    find_closest(freqs_prestim, multiple_50);
[~, multiple_50_ids_pststim] = ...
    find_closest(freqs_pststim, multiple_50);

% Replace 50Hz+-~3Hz multiple with nan 
% freq step =  0.6Hz for pre-stim
% freq step = 0.15Hz for post-stim
for freq_id = -20:20
    if freq_id > -6 && freq_id < 6
        pow_spect_prestim(multiple_50_ids_prestim+freq_id, :) = nan;
    end %freq_id > -6 && freq_id << 6
    pow_spect_pststim(multiple_50_ids_pststim+freq_id, :)= nan;        
end % freq_id = -10:10

    
% mean/std across trials
mean_prestim = nanmean(pow_spect_prestim, 2); 
mean_pststim = nanmean(pow_spect_pststim, 2);
std_prestim = nanstd(pow_spect_prestim, 0, 2);
std_pststim = nanstd(pow_spect_pststim, 0, 2);
stderr_prestim = std_prestim / sqrt(size(pow_spect_prestim, 2));
stderr_pststim = std_pststim / sqrt(size(pow_spect_pststim, 2));

figure();clf
hold on
plot(freqs_pststim, mean_pststim);
plot(freqs_prestim, mean_prestim);
% For shading
%lower_prestim = mean_prestim - stderr_prestim;
%lower_pststim = mean_pststim - stderr_pststim;
%upper_prestim = mean_prestim + stderr_prestim;
%upper_pststim = mean_pststim + stderr_pststim;
%h1 = fill([freqs_pststim fliplr(freqs_pststim)], ...
%          [upper_pststim fliplr(lower_pststim)], 'r');
%h2 = fill([freqs_prestim fliplr(freqs_prestim)], ...
%           [upper_prestim fliplr(lower_prestim)], 'r');
%set(h1, 'FaceColor', col)
%set(h1, 'EdgeColor', col)
hold off
xlim_upper = min(Ny_freq_prestim, Ny_freq_pststim); 
xlim([0 xlim_upper]);
legend({'post-stimulus', 'pre-stimulus'});
xlabel('frequency [Hz]');
ylabel('log power [\muV^2s]')

set(gca, 'FontSize', 10);
set(gca, 'FontName', 'Arial');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, 10, 10];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(f, img_fmt, 'FigRev3_prepost');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_prepost', '-r300');
end

% log narrow
xlim([0 500]);
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(f, img_fmt, 'FigRev3_prepost_narrow');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_prepost_narrow', '-r300');
end

% log Scale
xlim([0 xlim_upper]);
set(gca, 'XScale', 'log')
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(f, img_fmt, 'FigRev3_log');
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(f, img_fmt, 'FigRev3_log', '-r300');
end

%% bar plot difference within 50~150Hz multiples
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);
[~, ind_foi] = find_closest(freqs, foi);

n_order = 10;
start_ids = 50:50:50*n_order;
end_ids = 150:150:150*n_order;

for mult_id = 1:n_order
    start_id = start_ids(mult_id);
    end_id = end_ids(mult_id);
    % pre-stimulus
    [~, ids_prestim] = find_closest(freqs_prestim, [start_id, end_id]);
    pow_spect_pststim_tmp = nanmean(pow_spec_pstsim(start_id:end_id, :), 2);
    % post-stimulus
    [~, ids_pststim] = find_closest(freqs_pststim, [start_id, end_id]);
    mean_pststim = nanmean(pow_spect_pststim, 2);

end

function [powers, freq] = power_spect(data, fsample)
% compute power spectrum 
% Input: 
%   data: (time x trial)
%   fsample: sampling rate
% Output:
%   powers: power spectrum
%   freq: frequencies

% set chronux parameters if none provided
params.Fs = fsample;
params.tapers = [1 1];
params.pad = 1;

nTrials = size(data, 2);

% calculate by iteration
for trial = 1:nTrials
    % extract
    timetrial = squeeze(data(:, trial));
        
    % calculate
    [p_part,f_part] = mtspectrumc(timetrial, params);
        
    % store
    powers(:, trial) = 10*log10(p_part');  %#ok<AGROW> don't know how many freqs

    if trial==1
        freq = f_part';
    end
        
end
end


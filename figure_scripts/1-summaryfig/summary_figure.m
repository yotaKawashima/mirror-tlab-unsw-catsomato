%% Set overall variables

% PREREQUISITES:
% Add the git repo to path.

% directory where the drive is mounted. 
disk_path = '/media/rannee/UNSW_Cat_Somatos/';

% set figure path to be this temporary directory that is not on the main
% script or figure paths (this is a directory where the old and broken
% scripts are kept).
fig_path = [];

% highest level directory of data
data_path = [disk_path 'data/'];

% directory with the raw files
raw_path = [data_path 'raw/'];

% path to third party toolboxes
tool_path = [disk_path 'scripts/thirdparty_toolboxes/'];
chron_dir = [tool_path 'chronux/'];

%% Figure 1, part a
% Image of the cat brain. not implemented in matlab.

% Figure 1 common variables
fig1_catname = 'C20110808_R03';
fig1_area = 'S1';
fig1_rawfile = 'Cat20110808_UtahRight-3_LFPStimSegs.mat';

fsize = 16; % font size

%% Figure 1, part c
% a graph of the time domain 

% choose the time range
timeran = [-0.2, 0.5];
chan = 256;

% list the data directory.
data_dir = [data_path 'included_datasets/' fig1_catname '/epoched_rsampsl_biprref'];

% want the max stimulation condition
fnames = dir(fullfile(data_dir, [fig1_catname, '*', fig1_area, '*.mat'])); 
% find all files for this cat
nConds = length(fnames);
load(fullfile(data_dir, fnames(nConds).name)) % load the last one

% need to reset the 0 timepoint. To do this, need to load (part of) the raw. 
load([raw_path fig1_rawfile], 'stimTime')

% originally, the time is set to 0 when the ramp starts. so to set 0 where
% the vibration starts, need to subtract rampup and presine
offset = stimTime.rampup + stimTime.presine;
offset_time = data.time{1} - offset;

% find the domain of the plot
[~, minind] = find_closest(offset_time, timeran(1));
[~, maxind] = find_closest(offset_time, timeran(2));

% find mean and standard deviation
plot_mean = mean(data.trial(chan, minind:maxind, :)*10, 3);
plot_std = std(data.trial(chan, minind:maxind, :)*10, 0, 3);

plot_time = offset_time(minind:maxind);


% finally, time to plot.
col = [175,238,238]/256;

figure(13)
hold on
h=fill([plot_time fliplr(plot_time)], [plot_mean+plot_std fliplr(plot_mean-plot_std)], 'r');
set(h, 'FaceColor', col)
set(h, 'EdgeColor', col)
plot(plot_time, plot_mean, 'LineWidth', 2)
hold off
xlabel('time (sec)', 'FontSize', fsize)
ylabel('amplitude (\muV)', 'FontSize', fsize)
set(gca, 'FontSize', fsize)

% and print
print(13, '-depsc', [fig_path 'figure1c_20110808R03chan256.eps'])

%% Figure 1, part b
% A graph of the paradigm. Not implemented in matlab.

%% Figure 1, part d
% max stim condition, same time, frequency
% we don't actually have a dataset that is the full time period as a
% frequency so extract the data above

data_singchan = data;
data_singchan.trial = data.trial(chan, :, :);
data_singchan.custom.nsignals = 1;

% add chronux to path for power calculation
addpath(genpath(chron_dir))

data_out = chronux_pwrspec_v2(data_singchan);

% remove chronux as we don't need it any more
rmpath(genpath(chron_dir))

[~,fend] = find_closest(data_out.freq{1}, 250);

mean_freq = mean(data_out.trial(:, 1:fend, :), 3);
freq_data = data_out.freq{1}(1:fend);
figure(14)
plot(freq_data, mean_freq)
axis tight
xlabel('frequency (Hz)', 'FontSize', fsize)
ylabel('power (10log_{10}(\muV^s/s^2)', 'FontSize', fsize)
set(gca, 'FontSize', fsize)

print(14, '-depsc', [fig_path 'figure1d_20110808R03chan256.eps'])
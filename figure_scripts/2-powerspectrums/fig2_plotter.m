% Figure 2 plotter.

%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Setup
% file path
% dir_sec = fullfile(data_path, 'collated_data', 'anova_props');%'anoved_rsampsl_biprref_cmtspwr'
data_dir = fullfile(data_path, 'included_datasets');
save_dir = fullfile(data_path, 'collated_data', 'fig2_data');

% output image type
imgtype = '-depsc';

% Current date and time for output image
dt = datestr(now, 'yyyymmdd');

%% Options
max_only = true;
fmax = 250;
reload = false;

%% A: log power.

[data, ~] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr', max_only, reload, false, fmax);

figure(1); clf
plot(data.freq{1}, nanmean(data.trial, 1))

%% B/D: SNR

[data, ttest_data] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr', max_only, reload, true, fmax);

figure(2); clf
plot(data.freq{1}, nanmean(data.trial, 1))

figure(3); clf
plot(ttest_data.freq{1}, ttest_data.tstats)

%% C/E: EP

[data, ttest_data] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr', max_only, reload, true, fmax);

figure(4); clf
plot(data.freq{1}, nanmean(data.trial, 1))

figure(5); clf
plot(ttest_data.freq{1}, ttest_data.tstats)

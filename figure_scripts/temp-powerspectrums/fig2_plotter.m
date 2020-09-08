% Figure 2 plotter.
% Plot 5 values:
% 1. mean power across all trials
% 2. log SNR
% 3. t-statistic when using log SNR
% 4. evoked power
% 5. t-statistic when using evoked power


%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Setup
% file path
% dir_sec = fullfile(data_path, 'collated_data', 'anova_props');%'anoved_rsampsl_biprref_cmtspwr'
data_dir = fullfile(data_path, 'included_datasets');
save_dir = fullfile(data_path, 'collated_data', 'fig2_data');
% mkdir if not exist.
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% output image type
%imgtype = '-dpng';
imgtype = '-dsvg';

% Current date and time for output image
dt = datestr(now, 'yyyymmdd');

%% Options
max_only = true;
fmax = 250;
%reload = false;
reload = true;

if max_only
    m = 'max_';
else
    m = '';
end

%% A: log power.
% data mean across trials(repetitions) and channels 
[data, ~] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr', max_only, reload, false, fmax);

figure(1); clf
% data.trial dim = (channels*sessions) x frequencies
plot(data.freq{1}, nanmean(data.trial, 1)) % mean across channels&sessions
xlabel('frequency f [Hz]')
ylabel('log power [uV2S]')
print(gcf, imgtype, ['Fig2a_' m dt])

%% B/D: SNR

[data, ttest_data] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr', max_only, reload, true, fmax);

figure(2); clf
plot(data.freq{1}, nanmean(data.trial, 1)) 
xlabel('frequency f [Hz]')
ylabel('log SNR [-]')
print(gcf, imgtype, ['Fig2b_' m dt])

figure(3); clf
plot(ttest_data.freq{1}, ttest_data.tstats)
xlabel('frequency f [Hz]')
ylabel('t-statistic [-]')
print(gcf, imgtype, ['Fig2d_' m dt])

%% C/E: EP

[data, ttest_data] = f2p_loader(data_dir, save_dir, 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr', max_only, reload, true, fmax);

figure(4); clf
plot(data.freq{1}, nanmean(data.trial, 1)) % mean across channels & sessions
xlabel('frequency f [-]')
ylabel('log evoked power [-]')
print(gcf, imgtype, ['Fig2c_' m dt])

figure(5); clf
plot(ttest_data.freq{1}, ttest_data.tstats)
xlabel('frequency f [-]')
ylabel('t-statistic [-]')
print(gcf, imgtype, ['Fig2e_' m dt])

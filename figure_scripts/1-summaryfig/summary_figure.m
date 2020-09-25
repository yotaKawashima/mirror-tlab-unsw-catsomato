% Description
% Plot two values:
% 1. time course of bipolar-rereferenced voltage (mean across trials).
%    Bioplar channel 156 in between 28 electrode and 38 electrode.
%       You can check this by showing data.label(channelid).
%    Time from 0.2s to 0.5s. (time from stimulus onset.)
%    
% 2. power spectrum
%    Plot power spectrum at the channel.
%    
% 3. log SNR
%    Plot log SNR 
%
% 4. Vibration Evoked Log Power (VELogP)
%    Plot VELogP
% - note - 
% Inlcude fieldtrip-20140612 dir in path will give an error in
% mtspectrumc.m function.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../path_setup.m'))

[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

% image format
%img_fmt = '-dtiff';
img_fmt = '-dpdf';

% Figure 1 common variables
fig1_catname = 'C20110808_R03';
fig1_area = 'S1';

fsize = 10; % font size
ftype = 'Arial'; % font type
x_width = 6.5; % fig width
y_width = 4.5; % fig height

% Frequencies of interest (Just for visualisatoin)
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

%% Figure 1, part a
% Image of the cat brain. not implemented in matlab.

%% Figure 1, part b
% A graph of the paradigm. Not implemented in matlab.

%% Figure 1, part c
% a graph of the time domain (bipolar-rereferenced)
% Subtract baseline from each trial. Here we define a baseline per trial.
% i.e.
% baseline(trial) = mean(voltage(-0.2s:-0.1s,trial))
% Finally, we show vol(trial) = vol(:,trial) - baseline(trial)

% choose the time range
timeran = [-0.5, 2.5]; % from stimulus onset
timebaseline = [-0.5, -0.1]; % Baseline period
chan = 256; % Note that the first 100 channels are not bipolar channels. 

% list the data directory.
data_dir = [data_path 'included_datasets/' fig1_catname '/epoched_rsampsl_biprref'];

% want the max stimulation condition
fnames = dir(fullfile(data_dir, [fig1_catname, '*', fig1_area, '*.mat'])); 
nConds = length(fnames);
% load the last one (max stimulation condition)
load(fullfile(data_dir, fnames(nConds).name)) 

% find the domain of the plot
[~, minind] = find_closest(data.time{1}, timeran(1));
[~, maxind] = find_closest(data.time{1}, timeran(2));
% find the baseline period
[~, minindbaseline] = find_closest(data.time{1}, timebaseline(1));
[~, maxindbaseline] = find_closest(data.time{1}, timebaseline(2));

% Subtract baseline from voltage per trial
baselines = squeeze(mean(data.trial(chan, minindbaseline:maxindbaseline,:), 2));
bipolarvol = squeeze(data.trial(chan,minind:maxind,:));
subtvol = bipolarvol - baselines';

% Convert unit from micro V to milli V
subtvol = subtvol / 100;

% find mean and standard deviation
% unit of this 
%plot_mean = mean(data.trial(chan, minind:maxind, :)*10, 3);
%plot_std = std(data.trial(chan, minind:maxind, :)*10, 0, 3);
plot_mean = mean(subtvol, 2);
plot_std = std(subtvol, 0, 2); % unbiased

% Correct dim of matrix
plot_mean = plot_mean';
plot_std = plot_std'; 

plot_time = data.time{1}(minind:maxind);

% finally, time to plot.
%col = [175,238,238]/256; %Light blue
col = [192,192,192]/256;

figure(13)
hold on
h=fill([plot_time fliplr(plot_time)], [plot_mean+plot_std fliplr(plot_mean-plot_std)], 'r');
set(h, 'FaceColor', col)
set(h, 'EdgeColor', col)
%errorbar(plot_time, plot_mean, plot_std, 'k', 'LineWidth', 1)
plot(plot_time, plot_mean, 'k', 'LineWidth', 0.5)
hold off
xlim(timeran);
xlabel('time [s]', 'FontSize', fsize);
ylabel('amplitude [mV]', 'FontSize', fsize);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(13, img_fmt, [fig_path 'figure1_c']);
elseif img_fmt == "-dtiff"
    print(13, img_fmt, [fig_path 'figure1_c'], '-r300');
end

% Zoom in
%set(gca, 'Xlim', [0.5, 0.51]);
%print(13, '-dpng', [fig_path 'ZoomInFig1c_20110808R03chan256_' dt])

%% Figure 1, part d
% A graph of powerspectrum 
% max stim condition, same time, frequency
% Here, we show powerspectrum estimated from data from 0.5s to 2.5s after
% stimulus onset.

% Load powerspectrum data 
% list the data directory.
data_dir = [data_path 'included_datasets/' fig1_catname '/epoched_rsampsl_biprref_evkresp_cmtspwr'];

% Compure mean and std across trials for this channel data. 
[mean_logpower, std_logpower, freq_data] = mean_std_data(data_dir, fig1_catname, chan);

% Plot
figure(14)
hold on
h2=fill([freq_data fliplr(freq_data)], [mean_logpower+std_logpower fliplr(mean_logpower-std_logpower)], 'r');
set(h2, 'FaceColor', col)
set(h2, 'EdgeColor', col)
%errorbar(freq_data, mean_logpower, std_logpower, 'k', 'LineWidth', 1)
plot(freq_data, mean_logpower, 'k', 'LineWidth', 0.5)
axis tight 


ylim([0,inf])
xlabel('frequency [Hz]', 'FontSize', fsize)
ylabel('log power [\muV^2s]', 'FontSize', fsize)
% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [0, max(mean_logpower)*0.05];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 1)
xax = gca;
set(xax.XAxis, 'TickDir', 'out');
hold off
set(gca, 'XLim', [0,250]);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(14, img_fmt, [fig_path 'figure1_d']);
elseif img_fmt == "-dtiff"
    print(14, img_fmt, [fig_path 'figure1_d'], '-r300');
end


% Zoom in
%set(gca, 'Xlim', [100,110]);
%set(gca, 'Ylim', [0,30]);
%print(14, '-dpng', [fig_path 'ZoomInFig1d_20110808R03chan256_' dt])

%% Figure 1, part e
% A graph of log SNR
% max stim condition, same time, frequency
% Here, we show powerspectrum estimated from data from 0.5s to 2.5s after
% stimulus onset.

% Load powerspectrum data 
% list the data directory.
data_dir = [data_path 'included_datasets/' fig1_catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'];

% Compure mean and std across trials for this channel data. 
[mean_logsnr, std_logsnr, freq_data] = mean_std_data(data_dir, fig1_catname, chan);

% Plot 
figure(15)
hold on
h2=fill([freq_data fliplr(freq_data)], [mean_logsnr+std_logsnr fliplr(mean_logsnr-std_logsnr)], 'r');
set(h2, 'FaceColor', col);
set(h2, 'EdgeColor', col);
%errorbar(freq_data, mean_logsnr, std_logsnr, 'k', 'LineWidth', 1)
plot(freq_data, mean_logsnr, 'k', 'LineWidth', 0.5);
axis tight
%ylim([0,inf])
xlabel('frequency [Hz]', 'FontSize', fsize)
ylabel('log SNR [dB]', 'FontSize', fsize)
% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(1)+(ylims(2)-ylims(1))*0.05)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 1)
xax = gca;
set(xax.XAxis, 'TickDir', 'out');
hold off
set(gca, 'XLim', [0,250]);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(15, img_fmt, [fig_path 'figure1_e']);
elseif img_fmt == "-dtiff"
    print(15, img_fmt, [fig_path 'figure1_e'], '-r300');
end

% Zoom in
%set(gca, 'Xlim', [100,110]);
%print(15, '-dpng', [fig_path 'ZoomInFig1e_20110808R03chan256_' dt])

%% Figure 1, part f
% A graph of VELogP
% max stim condition, same time, frequency
% Here, we show powerspectrum estimated from data from 0.5s to 2.5s after
% stimulus onset.

% list the data directory.
data_dir = [data_path 'included_datasets/' fig1_catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr'];

% Compure mean and std across trials for this channel data. 
[mean_logvep, std_logvep, freq_data] = mean_std_data(data_dir, fig1_catname, chan);

% Plot 
figure(16)
hold on
h2=fill([freq_data fliplr(freq_data)], [mean_logvep+std_logvep fliplr(mean_logvep-std_logvep)], 'r');
set(h2, 'FaceColor', col)
set(h2, 'EdgeColor', col)
%errorbar(freq_data, mean_logvep, std_logvep, 'k', 'LineWidth', 1)
plot(freq_data, mean_logvep, 'k', 'LineWidth', 0.5)
axis tight
%ylim([0,inf])
xlabel('frequency [Hz]', 'FontSize', fsize)
ylabel('VELogP [dB]', 'FontSize', fsize)
% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(1)+(ylims(2)-ylims(1))*0.05)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 1)
xax = gca;
set(xax.XAxis, 'TickDir', 'out');
hold off
set(gca, 'XLim', [0,250]);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc"  || img_fmt == "-dpdf"   
    print(16, img_fmt, [fig_path 'figure1_f']);
elseif img_fmt == "-dtiff"
    print(16, img_fmt, [fig_path 'figure1_f'], '-r300');
end

% Zoom in
%set(gca, 'Xlim', [100,110]);
%print(16, '-dpng', [fig_path 'ZoomInFig1f_20110808R03chan256_' dt])


%% Functions
function [mean_data, std_data, freq_data] = mean_std_data(data_dir, fig1_catname, chan)
    % Load S1 data only max condition. Then compute mean and std across trials.
    % Input:    data_dir = directory name
    %           fig1_catname = cat name 
    %           chan = channels we plot
    % Output:   mean_data = mean across trials
    %           std_data = std across trials
    %           freq_data = frequencies (0-250Hz)

    % want the max stimulation condition
    fnames = dir(fullfile(data_dir, [fig1_catname, '*S1*.mat'])); 
    nConds = length(fnames);
    % load the last one (max stimulation condition)
    load(fullfile(data_dir, fnames(nConds).name))

    % Extract data from the channels (dim = frequencies x trials)
    trial_data = squeeze(data.trial(chan,:,:)); 

    % Get index of 250Hz
    [~,fend] = find_closest(data.freq{1}, 250);

    % Take mean of log VEP across trials.
    mean_data = mean(trial_data(1:fend, :), 2);
    std_data = std(trial_data(1:fend, :), 0, 2); % unbiased
    % Correct dim
    mean_data = mean_data';
    std_data = std_data';
    % Frequency data 
    freq_data = data.freq{1}(1:fend);
    freq_data = freq_data';
end

%% Previous versions
%% Figure 1, part d
% Rennee's version
% She included full time period here. i.e. 3s for this. 
% In the main script, we include only 2s data from 0.5s to 2.5s after
% stimulus onset.

%{
% A graph of powerspectrum 
% max stim condition, same time, frequency
% we don't actually have a dataset that is the full time period as a
% frequency so extract the data above. (We )

data_singchan = data;
data_singchan.trial = data.trial(chan, :, :);
data_singchan.custom.nsignals = 1;

% add chronux to path for power calculation
addpath(genpath(chron_dir))

% Compute power spectrum 
% Version 2 uses mtspectrumc's ability to compute a number of trials
% together for one channel. (Compute power for each trial.)
% This function returns 10log10(power).
data_out = chronux_pwrspec_v2(data_singchan);

% remove chronux as we don't need it any more
rmpath(genpath(chron_dir))

[~,fend] = find_closest(data_out.freq{1}, 250);

% Take mean of power across trials.
mean_freq = mean(data_out.trial(:, 1:fend, :), 3);
std_freq = std(data_out.trial(:, 1:fend, :), 0, 3); % unbiased
freq_data = data_out.freq{1}(1:fend);
figure(14)
hold on
h2=fill([freq_data' fliplr(freq_data')], [mean_freq+std_freq fliplr(mean_freq-std_freq)], 'r');
set(h2, 'FaceColor', col)
set(h2, 'EdgeColor', col)
plot(freq_data, mean_freq, 'LineWidth', 0.1)
hold off
axis tight
ylim([0,inf])
xlabel('frequency [Hz]', 'FontSize', fsize)
ylabel('log power [\muV^2s]', 'FontSize', fsize)
set(gca, 'FontSize', fsize)

print(14, img_fmt, [fig_path 'Fig1d_20110808R03chan256_' dt])
%}
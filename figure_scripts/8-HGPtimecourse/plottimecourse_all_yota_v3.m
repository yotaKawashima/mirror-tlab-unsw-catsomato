%% Plot time course
% Shorten timewindow from 2s to 0.5s for time course.
% This leads to precision trade-off between freq and time. 
% Half Band Width (HBW) increases 0.5Hz to 2Hz. (We use one taper.)
% HBW = (k + 1)/2T
% Take top 10% channels of results of ANOVA performed on logP (or VELogP).
% Mean f-stats across frequencies (from -50Hz to 150Hz) are used for
% the channel selection. (Almost the same results with the other methods.)
% This is used for the paper.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
%datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

% select datatype
datatype_cmtsgrm = 'epoched_rsampsl_biprref_cmtsgrm';
datatype_signif = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

%img_fmt = '-dtiff';
%img_fmt = '-dpng';
%img_fmt = '-depsc';
img_fmt = '-dpdf';

fsize = 10; % font size
ftype = 'Arial'; % font type
x_width = 12; % fig width
y_width = 7; % fig height

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));

hbw = 2; % half bandwidth

%% Load ANOVA results (f-stats) from VELogP. 
a = 2; % S2
area = ['S' num2str(a)]; %area 
    
f2main = [];
sessions = [];
bp_channel_id = [];

% load f-stats from all conditions and store it to allf.
for c = 1:numel(cat_names)
    % find file names
    c_path =  fullfile(data_dir, cat_names{c}, datatype_signif);
    loadname = dir(fullfile(c_path, ['*', area, '*.mat']));
    load(fullfile(c_path, loadname.name))
    % check if unipolar id included in srn data. If included,
    % unipolar should be removed. 
    if strcmp(metavars.chanlabel{1}(1:3), 'raw')
        % the first bipolar channel
        bipolarchid = 1 + prod(metavars.custom.spatialconfig);
        % # of bipolar channels
        nChan = metavars.custom.nsignals - bipolarchid + 1; 
    else
        bipolarchid = 1;
        nChan = metavars.custom.nsignals;
    end

    % loaded data: fstats = (channels, freqs, type of f-stats=3)
    % Type of f-stats consists of 
    % (F2 main effect, F1 main effect, interaction).

    % Get only f2 main effect at HGB from bipoolar channels
    f_tmp = squeeze(fstats(bipolarchid:end, :, 1));
    [~, f50_ind_fstats] = find_closest(metavars.freq{1, 1}, 50);
    [~, f150_ind_fstats] = find_closest(metavars.freq{1, 1}, 150);
    f_tmp = f_tmp(:, f50_ind_fstats:f150_ind_fstats);
    f2main = [f2main; f_tmp]; % (bipolar channels*sessions) x freq
    sessions = [sessions; ones(nChan, 1)*c];
    bp_ch_tmp = 1:nChan;
    bp_channel_id = [bp_channel_id; bp_ch_tmp'];        
end


%% Get top 10% channels (sessions x bipolar channels) based on VELogP
% Firstly, compute mean f-stats across frequencies from 50Hz to 150Hz.
% Then, take top 10% of channels based on the mean f-stats.

% Mean f stats across frequncies.
mean_fstats = mean(f2main, 2);

% Take top 10% channels based on the count.
top10_threshold = prctile(mean_fstats, 90, 1);
f2main_top10 = find(mean_fstats > top10_threshold);
sessions_top10 = sessions(f2main_top10);
bp_channel_id_top10 = bp_channel_id(f2main_top10);


%% Load power spectrum data.
n_timecourses = [];
trial_flag = [];
timeaxis_cell = cell(numel(cat_names), 1);

for c = 1:numel(cat_names)
    c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, datatype_cmtsgrm);
    loadname = dir(fullfile(c_path, ['*', area,'*.mat']));
    if isempty(loadname)
        % add the chronux toolbox
        addpath(genpath(chron_dir))

        % call function 
        % Input bipolar data path and recompute power spectrum with a 
        % shorter window.
        % This function recompute power spectrum for each session.
        % Note "timecourse" function always outputs only bipolar channels.
        % i.e. uni-poalr channels are not included.
        bipolar_dir = fullfile(data_path, 'included_datasets', ...
                               cat_names{c}, datatype_cmtsgrm(1:end-8), '/');
        timecourse(bipolar_dir, cat_names{c}, area)

        % remove the chronux toolbox
        rmpath(genpath(chron_dir))

        loadname = dir(fullfile(c_path, ['*', area, '*.mat']));
    end

    % Get max condition only
    loadname = loadname(end); 
    % Load data and store it.
    % Loaded data has bipolar channels x freqs x trials x time
    load(fullfile(c_path, loadname.name))

    % Get bipolar channels to be considered in this session.
    session_now = (sessions_top10 == c);
    iChans = bp_channel_id_top10(session_now);

    % Extract data from -0.5s to 2.5s. 
    % data.custom.time - time course of bipolar re-ref data. 
    % data.freq_t - time courese of power spectrum. 
    timeaxis = data.freq_t+data.custom.time(1); 
    [~, t1] = find_closest(timeaxis, -0.5);
    [~, t2] = find_closest(timeaxis,  2.5);
    timeaxis_cell{c,1} = timeaxis(t1:t2);

    % Extract data from -0.5s to 2.5s.
    data_timeaxis = data.trial(:, :, :, t1:t2);

    if ~isempty(iChans) % Store data if there are significant bipolar channels.
        data_ch = data_timeaxis(iChans, :, :, :);

        % Compute baseline per trial (-0.5s to -0.25s)
        % take mean across 10 log power -0.5s to -0.25s for each trial.        
        [~, bt1] = find_closest(timeaxis, -0.5);
        [~, bt2] = find_closest(timeaxis, -0.25);
        baseline_per_trial = mean(data_ch(:,:,:,bt1:bt2), 4);
        normalised_timecourse = data_ch - baseline_per_trial;

        % Add nan value trials for concatenation time course data
        % across (bipolar channels*sessions).
        % Some sessions have 10 trials and the others have 15 trials.
        if size(normalised_timecourse, 3) == 10
            nans = nan(size(normalised_timecourse, 1), ...
                       size(normalised_timecourse, 2), ...
                       5, size(normalised_timecourse, 4));
            normalised_timecourse = cat(3, normalised_timecourse, nans);    
            trial_flag = [trial_flag; ones(length(iChans), 1)];
        else
            trial_flag = [trial_flag; zeros(length(iChans), 1)];
        end 

        % n_timecourses = (bipolar channels*sessions)x freqs x trials x time
        n_timecourses = cat(1, n_timecourses, normalised_timecourse);

    end % if isempty(iChans) 

end    

%% Compute HGP
% Remark: time window = 0.5s, HBW = 2Hz
% Use from -0.5s to 2.5s from stimulus onset.

% Remove harmonics and intermodulation from 50Hz<f<150Hz
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

% Compute HGP 
% Get 50HZ<f<150Hz 
[mask, maskedf] = makefreqmask(data.freq{1}, foi, [50 150], hbw);
% Mean across freqs within HGB.
% i.e. channels x freqs x trials x time -> channels x trials x time
HGPs = squeeze(mean(n_timecourses(:, mask, :, :),2)); 
% Take mean across trials 
% i.e. channels x trials x time -> channels x time
mean_HGP_across_tri = squeeze(nanmean(HGPs, 2));
std_HGP_across_tri = squeeze(nanstd(HGPs, 0, 2));
% Compute sem across trials
sem_HGP_across_tri = ones(size(std_HGP_across_tri));
for i_ch = 1:size(std_HGP_across_tri, 1)
    if trial_flag(i_ch) == true % # of trials for this channel = 10
        sem_HGP_across_tri(i_ch, :) = std_HGP_across_tri(i_ch,:) / sqrt(10);
    else  % # of trials for this channel = 15
        sem_HGP_across_tri(i_ch, :) = std_HGP_across_tri(i_ch,:) / sqrt(15);            
    end
end

%% Compute f1 harmonics and intermodulation responses (each freq)
% Take responses at foi around HGP.
% find power at foi around HGP.
foi_around = foi(foi>40 & foi<160);
% Initialisation
% dim = fois x channels x times.
foi_resps_timecourse_across_t = zeros(length(foi_around), ...
                                      size(mean_HGP_across_tri, 1),...
                                      size(mean_HGP_across_tri, 2));
foi_resps_timecourse_std_across_t = zeros(length(foi_around), ...
                                      size(mean_HGP_across_tri, 1),...
                                      size(mean_HGP_across_tri, 2));
                                  
for i_foi = 1:length(foi_around)
    foi_now = foi_around(i_foi);
    [~, foi_ind] = find_closest(data.freq{1}, foi_now);
    % Take mean across trials for each channel for each foi
    % channels x freqs x trials x time -> channel x time
    foi_resps_timecourse_across_t(i_foi, :, :) = ...
        squeeze(nanmean(n_timecourses(:, foi_ind, :, :), 3)); 
    foi_resps_timecourse_std_across_t(i_foi, :, :) = ...
        squeeze(nanstd(n_timecourses(:, foi_ind, :, :), 0, 3)); 
end %i=1:length(foi_around)

%% Compute f1 harmonics and intermodulation responses (each type)
% Take mean across foi for each type. (harmonic or IM).
% type x channels x times
foi_resps_timecourse_type_mean_across_t = zeros(2, ...
                                                size(mean_HGP_across_tri, 1),...
                                                size(mean_HGP_across_tri, 2));
foi_resps_timecourse_type_std_across_t = zeros(2, ...
                                               size(mean_HGP_across_tri, 1),...
                                               size(mean_HGP_across_tri, 2));
foi_resps_timecourse_type_sem_across_t = zeros(2, ...
                                               size(mean_HGP_across_tri, 1),...
                                               size(mean_HGP_across_tri, 2));
for i_type = 1:2
    switch i_type 
        case 1 % f1 harmonics 
            foi_now = foi_f1_and_harm;
        case 2 % IM
            foi_now = foi_inter;
    end
    foi_now = foi_now(foi_now>40 & foi_now<160);
    [~, foi_now_inds] = find_closest(data.freq{1}, foi_now);
    % Compute mean across frequencies (Similar to HGP definition)
    % channels x freqs x trials x time -> channels x trials x time
    n_timecourses_now = squeeze(mean(n_timecourses(:, foi_now_inds, :, :), 2));

    % Mean and std across trials
    % channels x trials x time -> channels x time
    foi_resps_timecourse_type_mean_across_t(i_type, :, :) = squeeze(nanmean(n_timecourses_now, 2));
    foi_resps_timecourse_type_std_across_t(i_type, :, :) = squeeze(nanstd(n_timecourses_now, 0, 2));
    % SEM across trials
    for i_ch = 1:size(foi_resps_timecourse_type_std_across_t, 2)
        if trial_flag(i_ch) == true % # of trials for this channel = 10
            foi_resps_timecourse_type_sem_across_t(i_type, i_ch, :) = ...
                    foi_resps_timecourse_type_std_across_t(i_type, i_ch, :) / sqrt(10);
        else  % # of trials for this channel = 15
            foi_resps_timecourse_type_sem_across_t(i_type, i_ch, :) = ... 
                    foi_resps_timecourse_type_std_across_t(i_type, i_ch, :) / sqrt(15);            
        end
    end

end % i_type = 1:2

%% Show # of channels showing F2 main effect significance.
fprintf('# of channels: %i\n', size(f2main_top10, 1));

%% Plot time course HGP and each foi separately (Mean across channels*trials)
% Remark: time window = 0.5s, HBW = 2Hz
% Note time axis are slightly different across sessions due to the
% difference of timestamp. Here, I ignored the differences.

figure();clf
hold on

% Plot tagged power around HGB.
for i_foi=1:length(foi_around)
    foi_now = foi_around(i_foi);
    if sum(ismember(foi_f1_and_harm, foi_now)) % F1 harmonics
        color_ = [255, 47, 64]/256;
    elseif sum(ismember(foi_inter, foi_now)) % Intermodulation     
        color_ = [0, 82, 255]/256;
    end % sum(ismember(foi_f1_and_harm, foi_now)) 

    switch ceil(i_foi/2)
        case 1; linestyle = '-'; linewidth = 1;
        case 2; linestyle = '--'; linewidth = 1;
        case 3; linestyle = ':'; linewidth = 1;
        case 4; linestyle = '-.'; linewidth = 1;
        case 5; linestyle = '-'; linewidth = 2;
    end % switch ceil(i/3)
    % Mean across channels 
    mean_across_c = squeeze(mean(foi_resps_timecourse_across_t(i_foi,:,:),2));
    plot(timeaxis_cell{1,1}, mean_across_c, 'Color', color_, ...
        'linewidth', linewidth, 'linestyle', linestyle);            
end % i=1:length(foi_around)

% Plot HGP time course
mean_HGP = squeeze(mean(mean_HGP_across_tri, 1)); % Mean across channels
plot(timeaxis_cell{1,1}, mean_HGP, 'color', 'k', 'linewidth', 2);   


%Show y=0 line
hline = yline(0, 'color', 'k', 'LineWidth', 0.5);
% Remove legend for these lines
hline_annotation = get(hline, 'Annotation');
set(get(hline_annotation,'LegendInformation'), 'IconDisplayStyle','off');
hold off

% Setting 
xlim([-0.5, 2.5]);
ylim([-4, 10]);
legend_tagged = cellstr(strcat(num2str(foi_around'), ' Hz'))';
legend([legend_tagged, {'HGP'}], 'Location',  'eastoutside');
xlabel('time [s]');
ylabel('normalised power [dB]');
set(gca,'position',[0.1 0.2, 0.5, 0.7]);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');    
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, [fig_path 'figure8_a']);
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, [fig_path 'figure8_a'], '-r300');
end

%% Plot timecourse HGP, harmonics, and IM (mean and SEM across channels*trials)
figure();clf

hold on
% Plot tagged power around HGB.
p = gobjects(3);
for i_type = 1:2

    if i_type ==1  % F1 harmonics
        color_ = [255, 47, 64]/256;
    elseif i_type ==2 % Intermodulation     
        color_ = [0, 82, 255]/256;
    end % if i_type == 1
    % Mean across channels
    mean_now = squeeze(mean(foi_resps_timecourse_type_mean_across_t(i_type, :, :), 2))';
    % mean sem across channels
    sem_now = squeeze(mean(foi_resps_timecourse_type_sem_across_t(i_type, :, :), 2))';

    p(i_type)= plot(timeaxis_cell{1,1}, mean_now,...
                    'Color', color_, 'linewidth', 1.5);
    h = fill([timeaxis_cell{1,1} fliplr(timeaxis_cell{1,1})],...
             [mean_now+sem_now fliplr(mean_now-sem_now)], 'r');
    set(h, 'FaceColor', color_);
    set(h, 'EdgeColor', color_);
    % Changing alpha causes loading eps file error.
    set(h, 'FaceAlpha', 0.3); 
    set(h, 'EdgeAlpha', 0.3);

end % for i_type = 1:2

% Mean HGP and its sem across channels
mean_HGP = squeeze(mean(mean_HGP_across_tri, 1)); % Mean across channels
mean_sem_HGP = squeeze(mean(sem_HGP_across_tri, 1)); 
%uistack(h1, 'bottom');
p(3) = plot(timeaxis_cell{1,1}, mean_HGP, 'color', 'k', 'linewidth', 1.5);
h = fill([timeaxis_cell{1,1} fliplr(timeaxis_cell{1,1})],...
         [mean_HGP+mean_sem_HGP fliplr(mean_HGP-mean_sem_HGP)], 'r');
set(h, 'FaceColor', [200, 200, 200]/256);
set(h, 'EdgeColor', [200, 200, 200]/256);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0.5);
%Show y=0 line
hline = yline(0, 'color', 'k', 'LineWidth', 0.5);
% Remove legend for these lines
hline_annotation = get(hline, 'Annotation');
set(get(hline_annotation,'LegendInformation'), 'IconDisplayStyle','off');
hold off

uistack(p(3), 'bottom');
uistack(p(2), 'bottom');
uistack(p(1), 'bottom');

% Setting 
xlim([-0.5, 2.5]);
ylim([-4, 10]);
legend({'F1 Harmonics', 'Intermodulations', 'HGP'}, 'Location', 'eastoutside');
xlabel('time [s]');
ylabel('normalised power [dB]');
set(gca,'position',[0.1, 0.2, 0.5, 0.7]);
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');    
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, [fig_path 'figure8_b']);
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, [fig_path 'figure8_b'], '-r300');
end

%print(gcf, img_fmt, ['HGPtimecourseline_V3_MEAN_SEM_' area '_' dt]);

%% Plot timecourse only some frequencies with std.
%{
figure();clf
hold on

% Plot tagged power around HGB.
foi_selected = [69, 85];

p = gobjects(3);
for i_type=1:length(foi_selected)
    %plot(timeaxis_cell{1,1}, foi_resps_timecourse(i,:));    

    foi_now = foi_selected(i_type);
    foi_ind = find(foi_around == foi_now);

    if sum(ismember(foi_f1_and_harm, foi_now)) % F1 harmonics
        color_ = [255, 47, 64]/256;
    elseif sum(ismember(foi_inter, foi_now)) % Intermodulation     
        color_ = [0, 82, 255]/256;
    end % sum(ismember(foi_f1_and_harm, foi_now)) 

    mean_now = foi_resps_timecourse(foi_ind, :);
    std_now = foi_resps_timecourse_std(foi_ind, :);

    p(i_type)= plot(timeaxis_cell{1,1}, mean_now,...
         'Color', color_, 'linewidth', 1.5,...
         'linestyle', '--');  

    h = fill([timeaxis_cell{1,1} fliplr(timeaxis_cell{1,1})],...
             [mean_now+std_now fliplr(mean_now-std_now)], 'r');

    set(h, 'FaceColor', color_);
    set(h, 'EdgeColor', color_);
    set(h, 'FaceAlpha', 0.1);
    set(h, 'EdgeAlpha', 0.1);

end % i_type=1:length(foi_around)

% Plot HGP time course
p(3) = plot(timeaxis_cell{1,1}, mean_HGP, 'color', 'k', 'linewidth', 3);

h = fill([timeaxis_cell{1,1} fliplr(timeaxis_cell{1,1})],...
         [mean_HGP+std_HGP fliplr(mean_HGP-std_HGP)], 'r');
set(h, 'FaceColor', [200, 200, 200]/256);
set(h, 'EdgeColor', [200, 200, 200]/256);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0.5);

hold off

uistack(p(3), 'bottom');
uistack(p(2), 'bottom');
uistack(p(1), 'bottom');

% Setting 
xlim([-0.5, 2.5]);

legend_tagged = cellstr(strcat(num2str(foi_selected'), ' Hz'))';
legend([legend_tagged, {'HGP'}]);
set(gca, 'fontsize', 16);
xlabel('time [s]');
ylabel('normalised power [dB]');
%}


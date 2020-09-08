%% Plot time course
% Shorten timewindow from 2s to 0.5s for time course.
% This leads to precision trade-off between freq and time. 
% Half Band Width (HBW) increases 0.5Hz to 2Hz. (We use one taper.)
% HBW = (k + 1)/2T
% Channels are relected based on AVNOVA results on HGP.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%% Set script specific variables
% select area
area = 'S2';

% select datatype
data_type = 'epoched_rsampsl_biprref_cmtsgrm';
data_type_signif = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp_adatain_adatout';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));

q = 0.05; % significance level for fdr
hbw = 2; % half bandwidth

%% Load ANOVA results for HGP. 
% load the p-values
ss_chans = cell(1, numel(cat_names));
s_chans = [];
for c = 1:numel(cat_names)
    c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, data_type_signif);
    loadname = dir(fullfile(c_path, '*S2*.mat'));
    load(fullfile(c_path, loadname.name))
    % check if unipolar id included in srn data. If included,
    % unipolar should be removed. 
    if strcmp(metavars.label{1}(1:3), 'raw')
        % the first bipolar channel
        bipolarchid = 1 + prod(metavars.custom.spatialconfig);
        % # of bipolar channels
        nChan = metavars.custom.nsignals - bipolarchid + 1; 
    else
        nChan = metavars.custom.nsignals;
        bipolarchid = 1;
    end
    % loaded data: pvals = (channels, nfreqs=1, type of f-stats=3)
    % Type of f-stats consists of 
    % (F2 main effect, F1 main effect, interaction).
   
    % Pick up only F2 main effect.
    s_tmp = pvals(bipolarchid:end, 1, 1); 
    s_chans = [s_chans; s_tmp]; % Pool all sessions for fdr.
    ss_chans{c} = s_tmp; % Store p-val for each session. 
end

%% fdr correction (pool all sessions and bipolar channels in S2.)
[pID, ~] = eeglab_fdr(s_chans, q, 'parametric');
% iChans = find(s_chans < pID);
% p_tmp = p_tmp(iChans, :, :, :);

%% Load power spectrum data.
n_timecourses = [];
timeaxis_cell = cell(numel(cat_names),1);
for c = 1:numel(cat_names)
    c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, data_type);
    loadname = dir(fullfile(c_path, ['*', area, '*.mat']));
    loadname = loadname(end); % get max condition only
    if isempty(loadname)
        % add the chronux toolbox
        addpath(genpath(chron_dir))
        
        % call function 
        % Input bipolar data path and recompute power spectrum with a 
        % shorter window.
        % This function recompute power spectrum for each session.
        % Note "timecourse" function always outputs only bipolar channels.
        % i.e. uni-poalr channels are not included.
        timecourse(c_path(1:end-8), cat_names{c}, area)
        
        % remove the chronux toolbox
        rmpath(genpath(chron_dir))
        
        loadname = dir(fullfile(c_path, ['*', area, '*.mat']));
        loadname = loadname(end); % get max condition only
    end
    % Load data and store it.
    % Loaded data has bipolar channels x freqs x trials x time
    load(fullfile(c_path, loadname.name))
    
    iChans = find(ss_chans{c} < pID); % Only significant channels
    
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
        
        % Permute into (bipolar channels*trials) x freqs x time
        for i=1:size(normalised_timecourse, 3)
            n_timecourses = [n_timecourses; 
                             squeeze(normalised_timecourse(:,:,i,:))];
        end
        
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
HGPs = squeeze(mean(n_timecourses(:, mask, :),2)); % Mean across freqs within HGB.
% Take mean across (bipolar channels*trials)
mean_HGP = squeeze(mean(HGPs, 1));
std_HGP = squeeze(std(HGPs,0,1));

% Take responses at foi around HGP.
% find power at foi around HGP.
foi_around = foi(foi>40 & foi<160);
foi_resps_foival = zeros(length(foi_around),1);
foi_resps_timecourse = zeros(length(foi_around),size(mean_HGP, 2));
foi_resps_timecourse_std = zeros(length(foi_around),size(mean_HGP, 2));
for i = 1:length(foi_around)
    foi_now = foi_around(i);
    [val, foi_ind] = find_closest(data.freq{1}, foi_now);
    foi_resps_foival(i,1) = val;
    % Take mean across (bipolar channels*trials)
    foi_resps_timecourse(i,:) = squeeze(mean(n_timecourses(:, foi_ind, :), 1))'; 
    foi_resps_timecourse_std(i,:) = squeeze(std(n_timecourses(:, foi_ind, :), 0, 1))'; 
end %i=1:length(foi_around)

%% Plot time course
% Remark: time window = 0.5s, HBW = 2Hz
% Note time axis are slightly different across sessions due to the
% difference of timestamp. Here, I ignored the differences.

figure();clf
hold on

% Plot tagged power around HGB.
for i=1:length(foi_around)
    %plot(timeaxis_cell{1,1}, foi_resps_timecourse(i,:));    
    foi_now = foi_around(i);
    if sum(ismember(foi_f1_and_harm, foi_now)) % F1 harmonics
        color_ = [255, 47, 64]/256;
    elseif sum(ismember(foi_inter, foi_now)) % Intermodulation     
        color_ = [0, 82, 255]/256;
    end % sum(ismember(foi_f1_and_harm, foi_now)) 

    switch ceil(i/2)
        case 1; linestyle = '-'; linewidth = 1.5;
        case 2; linestyle = '--'; linewidth = 1.5;
        case 3; linestyle = ':'; linewidth = 1.5;
        case 4; linestyle = '-.'; linewidth = 1.5;
        case 5; linestyle = '-'; linewidth = 3;
    end % switch ceil(i/3)
    plot(timeaxis_cell{1,1}, foi_resps_timecourse(i,:),...
         'Color', color_, 'linewidth', linewidth,...
         'linestyle', linestyle);             
end

% Plot HGP time course
plot(timeaxis_cell{1,1}, mean_HGP, 'color', 'k', 'linewidth', 3);

hold off

% Setting 
xlim([-0.5, 2.5]);
legend_tagged = cellstr(strcat(num2str(foi_around'), ' Hz'))';
legend([legend_tagged, {'HGP'}]);
set(gca, 'fontsize', 16);
xlabel('time [s]');
ylabel('normalised power [dB]');

% Print
print(gcf, '-dpng', ['HGPtimecourseline_V1_' area '_' dt]);

%% Show # of channels showing F2 main effect significance.
fprintf('# of F2 main effect significant channels: %i\n', sum(s_chans<pID));

%% Plot time course
% Shorten timewindow from 2s to 0.5s for time course.
% This leads to precision trade-off between freq and time. 
% Half Band Width (HBW) increases 0.5Hz to 2Hz. (We use one taper.)
% HBW = (k + 1)/2T
% Take top 10% channels of results of ANOVA performed on logP (or VELogP).
% Channels are selected based on how many times a channel can be top 10%
% across frequenices. (Almost the same results with the other methods.)

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
%datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

% select area
area = 'S2';

% select datatype
datatype_bipolar = 'epoched_rsampsl_biprref_cmtsgrm';
datatype_signif = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));

hbw = 2; % half bandwidth

%img_fmt = '-dpng';
img_fmt = '-depsc';
%% Load ANOVA results (f-stats) from VELogP. 
for a = 2:2
    area = ['S' num2str(a)]; %area 
    
    f2main = [];
    sessions = [];
    bp_channel_id = [];

    % load f-stats from all conditions and store it to allf.
    for c = 1:numel(cat_names)
        % find file names
        c_path =  fullfile(data_dir, cat_names{c}, datatype_signif);
        loadname = dir(fullfile(c_path, '*S2*.mat'));
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
    % Firstly, get channels giving f-stats (f2 main) above top 10% per
    % frequency within HGB (50Hz~150Hz). Secondly, check how many times 
    % each channel is counted in top 10% across frequencies. The more a 
    % channel is counted in, the more the channel is thought to contribute 
    % HGP responses. Finally, take top 10% channels based on the count.
    
    % Get top 10% channels per grequency
    fstats_prct10_per_freq = prctile(f2main, 90, 1); % per frequencies
    f2main_top10_per_freq = (f2main > fstats_prct10_per_freq);
    
    % Count how many time each channel is counted in as top 10%.
    count_across_freq = sum(f2main_top10_per_freq, 2);
    
    % Take top 10% channels based on the count.
    top10_channels = prctile(count_across_freq, 90, 1);
    f2main_top10 = find(count_across_freq > top10_channels);
    sessions_top10 = sessions(f2main_top10);
    bp_channel_id_top10 = bp_channel_id(f2main_top10);
       

    %% Load power spectrum data.
    n_timecourses = [];
    timeaxis_cell = cell(numel(cat_names), 1);
    
    for c = 1:numel(cat_names)
        c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, datatype_bipolar);
        loadname = dir(fullfile(c_path, ['*', area,'*.mat']));
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
            timecourse(c_path(1:end-7), cat_names{c}, area)

            % remove the chronux toolbox
            rmpath(genpath(chron_dir))

            loadname = dir(fullfile(c_path, ['*', area, '*.mat']));
            loadname = loadname(end); % get max condition only
        end
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
    legend([legend_tagged, {'HGP'}])
    set(gca, 'fontsize', 16);
    xlabel('time [s]');
    ylabel('normalised power [dB]');

    % Print
    print(gcf, '-dpng', ['HGPtimecourseline_V2_' area '_' dt]);

    %% Show # of channels showing F2 main effect significance.
    fprintf('# of channels: %i\n', size(f2main_top10, 1));

end

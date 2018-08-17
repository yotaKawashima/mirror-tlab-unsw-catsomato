%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Set script specific variables
% select area
S = 'S2';

% select datatype
data_type = 'epoched_rsampsl_biprref_cmtsgrm';
data_type_signif = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp_adatain_adatout';

% select condition
cond = 'F023A159_F200A016';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));
cat_names = cat_names(3);

% significance level
q = 0.05;
%% Load data
% load the p-values
s_chans = [];
for c = 1:numel(cat_names)
    c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, data_type_signif);
    loadname = dir(fullfile(c_path, '*S2*.mat'));
    
    load(fullfile(c_path, loadname.name))
    s_chans = [s_chans; permute(pvals(:, 65:end, 2), [2,3,1])]; % only get F2 channels
end

p_tmp = [];
for c = 1:numel(cat_names)
    c_path =  fullfile(data_path, 'included_datasets', cat_names{c}, data_type);
    loadname = dir(fullfile(c_path, ['*S2*' cond '*.mat']));
    if isempty(loadname) % TODO: not all cats have this data yet
        % add the chronux toolbox
        addpath(genpath(chron_dir))
        
        % call function
        timecourse(c_path(1:end-8), cat_names{c}, 'S2')
        
        % remove the chronux toolbox
        rmpath(genpath(chron_dir))
        
        loadname = dir(fullfile(c_path, ['*S2*' cond '*.mat']));
    end
    load(fullfile(c_path, loadname.name))
    p_tmp = [p_tmp; data.trial];
end

%% Get only the significant channels
[pID, ~] = eeglab_fdr(s_chans, q, 'parametric');
iChans = find(s_chans < pID);
p_tmp = p_tmp(iChans, :, :, :);

%% Plot colour images
% grab power data, mean across trials
% p_tmp = data.trial;
p_tmp = mean(p_tmp, 3); % channels x frequency x trials x timesteps
p_tmp = mean(p_tmp, 1);

% find frequency indices
[~, in] = find_closest(data.freq{1}, 100);
[~, i50] = find_closest(data.freq{1}, 50);
[~, i150] = find_closest(data.freq{1}, 150);
f_inds = i50+5:i150-5; % only concerned with the high gamma band



%% Plot line timecourse
% Remove harmonics and intermodulation
f1 = 23;
f2 = 200;
f_lo = 0;
f_hi = 250;
foi = sort(unique([f1:f1:f_hi f2:-f1:f_lo f2:f1:f_hi]));
[mask, maskedf] = makefreqmask(data.freq{1}, foi, [50.5 149.5], 0.5);

% find power outside of HGB
[mask2, ~] = makefreqmask(data.freq{1}, foi, [0.5 49.5], 0.5);
[mask3, ~] = makefreqmask(data.freq{1}, foi, [150.5 250], 0.5);

% find the baseline (from -0.7 to -0.2)
timeaxis = data.freq_t+data.custom.time(1);
[~, t1] = find_closest(timeaxis, -0.7);
[~, t2] = find_closest(timeaxis, -0.2);
p_tcf = permute(p_tmp, [4, 1, 2, 3]); % time x chan x freq
basefreq = mean(p_tcf(t1:t2, :, :), 1);
p_norm = p_tcf - repmat(basefreq, [size(p_tcf, 1),1,1]);

% plot
figure(3)
p = plot(timeaxis, mean(p_norm(:, :, mask),3), 'b', ...
    timeaxis, mean(p_norm(:, :, or(mask2,mask3)),3), 'r');
legend({'High gamma band', 'outside HGB'})
xlabel('time (s)')
ylabel('normalised power (dB)')

set(p(1), 'LineWidth', 3)
set(p(2), 'LineWidth', 2)

% TODO: re-enable printing
% print(gcf, '-dpng', ['HGPtimecourseline_' cat_name '_S' S '_Ch' num2str(chs(k), '%03i') '_' dt])

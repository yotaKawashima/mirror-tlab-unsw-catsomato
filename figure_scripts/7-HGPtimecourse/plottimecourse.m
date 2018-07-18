%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Set script specific variables
% select area
S = 'S2';

% select cat
cat_name = 'C20110808_R03';

% select condition
cond = 'F023A159_F200A016';

%% Load data

% get filename
data_type = 'epoched_rsampsl_biprref_cmtsgrm';
filename = [cat_name '_' S '_' cond '_' data_type '.mat'];
data_dir = fullfile(data_path, 'included_datasets', cat_name, data_type);
loadname = fullfile(data_dir, filename);

% find out if the data is already produced. if not, produce it.
if ~exist(loadname, 'file')
    % add the chronux toolbox
    addpath(genpath(chron_dir))
    
    % call function
    timecourse(data_dir(1:end-8), cat_naclme, S)
    
    % remove the chronux toolbox
    rmpath(genpath(chron_dir))
end

% load data
load(loadname)

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

%% Plot colour images
% grab power data, mean across trials
p_tmp = data.trial;
p_tmp = mean(p_tmp, 3); % channels x frequency x trials x timesteps

% find frequency indices
[~, in] = find_closest(data.freq{1}, 100);
[~, i50] = find_closest(data.freq{1}, 50);
[~, i150] = find_closest(data.freq{1}, 150);
f_inds = i50+5:i150-5; % only concerned with the high gamma band

% Find out which channels would be good
figure(1); shg; clf
imagesc(permute(mean(p_tmp(:, f_inds, :, :), 2), [1, 4, 2, 3]))

% Set channels to plot
if S == 'S1'
	chs = [174, 183, 235, 134];
else
    chs =102;%108;%[100, 125, 126, 96, 113];
end

% Plot for each channel
for k = 1:numel(chs)
    figure(2)
    imagesc(data.freq_t+data.custom.time(1), data.freq{1}(f_inds), permute(p_tmp(chs(k), f_inds, :, :), [2, 4, 1, 3])); 
    set(gca, 'ydir', 'normal'); 
    
    xlabel('time (s)')
    ylabel('f (Hz)')
    title(['C20110808_R08, ' S ' , channel ' num2str(chs(k))], 'interpret', 'none')
    
    colorbar
    
    shg
    
    % Save to file
    %print(gcf, '-dpng', ['TimeFreqFig_C20110808_R03_S' S '_Ch' num2str(chs(k), '%03i') '_' dt])
    %pause() % Wait for manual verification
end

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
plot(timeaxis, mean(p_norm(:, chs, mask),3), 'b', timeaxis, mean(p_norm(:, chs, or(mask2,mask3)),3), 'r')
legend({'High gamma band', 'outside HGB'})
xlabel('time (s)')
ylabel('normalised power (dB)')
print(gcf, '-dpng', ['HGPtimecourseline_' cat_name '_S' S '_Ch' num2str(chs(k), '%03i') '_' dt])
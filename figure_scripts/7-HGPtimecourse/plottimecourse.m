% select area
S = '2';

% load data
load(['C20110808_R03_S' S '_F023A159_F200A016_epoched_rsampsl_biprref_cmtsgrm.mat'])

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

% grab power data, mean across trials
p_tmp = data.trial;
p_tmp = mean(p_tmp, 3); % channels x frequency x trials x timesteps

% find frequency indices
[~, in] = find_closest(data.freq{1}, 100);
[~, i50] = find_closest(data.freq{1}, 50);
[~, i150] = find_closest(data.freq{1}, 150);
f_inds = i50+5:i150-5; % only concerned with the high gamma band

% Find out which channels would be good
figure(1); shg
imagesc(permute(mean(p_tmp(:, f_inds, :, :), 2), [1, 4, 2, 3]))

% Set channels to plot
if S == '1'
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
    title(['C20110808_R08, S' S ' , channel ' num2str(chs(k))], 'interpret', 'none')
    
    colorbar
    
    shg
    
    % Save to file
    %print(gcf, '-dpng', ['TimeFreqFig_C20110808_R03_S' S '_Ch' num2str(chs(k), '%03i') '_' dt])
    %pause() % Wait for manual verification
end


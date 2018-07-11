S = '2';
load(['C20110808_R03_S' S '_F023A159_F200A016_epoched_rsampsl_biprref_evkresp_cmtsgrm.mat'])
dt = datestr(now, 'yyyymmdd');


p_tmp = data.trial;
p_tmp = mean(p_tmp, 3);
[~, in] = find_closest(data.freq{1}, 100);
[~, i50] = find_closest(data.freq{1}, 50);
[~, i150] = find_closest(data.freq{1}, 150);

f_inds = i50+5:i150-5;

figure(1); shg
imagesc(permute(mean(p_tmp(:, f_inds, :, :), 2), [1, 4, 2, 3]))

if S == '1'
	chs = [174, 183, 235, 134];
else
    chs =[14, 34, 35, 94, 101, 102]+64;%108;%[100, 125, 126, 96, 113];
end

for k = 1:numel(chs)
    figure(2)
    imagesc(data.freq_t, data.freq{1}(f_inds), permute(p_tmp(chs(k), f_inds, :, :), [2, 4, 1, 3])); 
    set(gca, 'ydir', 'normal'); 
    
    xlabel('time (s)')
    ylabel('f (Hz)')
    title(['C20110808_R08, S' S ' , channel ' num2str(chs(k))], 'interpret', 'none')
    
    colorbar
    
    shg
    
    print(gcf, '-dpng', ['TimeFreqFig_C20110808_R03_S' S '_Ch' num2str(chs(k), '%03i') '_' dt])
    pause()
end


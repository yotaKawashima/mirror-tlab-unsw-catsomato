% concatenated f-statistic plot

% first define what data is being used. 
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';
dir_sec = [data_dir, data_type];

filenames = {'C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'};
% this version uses the 9 files that are also being analysed in the ANOVA
% results. 

nRec = numel(filenames);

options.load = false;
options.plotchannelscan = false;
options.plotchannel = {[]; []};
options.printtype = {'-dpng'};

%% run fstat_plot
for r = 1:nRec
    
    [mvars, f_stats] = fstat_plot(dir_sec, filenames{r}, options, [], []);
    
    if r == 1
        ms(1).m = mvars;
        fs = f_stats;
    else 
        ms(r).m = mvars;
        fs{1} = vertcat(fs{1},f_stats{1});
        fs{2} = vertcat(fs{2},f_stats{2});
    end
end

%% plot

area = 1;

av_fstat = mean(fs{area}, 1, 'omitnan'); % first dimension is channels
av_fstat = permute(av_fstat, [2, 3, 1]);

frange = [50 150];

stimf = [23, 200];

% remove harmonic frequencies
[~, keep_i] = no_resp_indices(stimf, frange, ms(1).m(1).freq{1});

length_freqs = keep_i(end)-keep_i(1)+1;
av_fstat_highgamma = NaN(length_freqs, size(av_fstat, 2));
av_fstat_highgamma(keep_i-keep_i(1)+1, :) = av_fstat(keep_i, :);

%figure
plot(ms(1).m(1).freq{1}(1:length_freqs), av_fstat_highgamma)
xlabel('frequency (Hz)')
ylabel('f-statistic')
legend({'200 Hz', '23 Hz', 'interaction'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')
title(['F-statistic against frequency in S' num2str(area) ', averaged across datasets'])
shg

print(gcf, '-dpng', ['concatenated_all_S' num2str(area) '_meanf-statplot'])
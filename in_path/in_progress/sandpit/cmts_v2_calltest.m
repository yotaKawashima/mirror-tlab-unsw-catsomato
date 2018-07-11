
filename = 'C20110808_R03';




data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/C20110808_R03/epoched_rsampsl_biprref_evkresp/';

load([data_dir 'C20110808_R03_TStim_S1_F023A159_F200A016_epoched_rsampsl_biprref_evkresp.mat'])

addpath(genpath('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/thirdparty_toolboxes/chronux'))

print_format = '-dpng';

taps = {[1, 1], [3, 5], [5, 9]};

for k = 1:numel(taps)
    params.Fs = data.fsample;
    params.tapers = taps{k};
    params.pad = 1;
    
    
    %data_out(k) = chronux_pwrspec_v2(data, params);
    
    f_ind = find_closest_ind(data_out(k).freq{1}, 300);
    
    figure(k)
    plot(data_out(k).freq{1}(1:f_ind), squeeze(mean(data_out(k).trial(256, 1:f_ind, :), 3)))
    print(gcf, print_format, ['tapers_test_' num2str(taps{k}(1)) '-' num2str(taps{k}(2)) '_power'])
    
    %snr_dat(k) = snr_1to3([], [], data_out(k));
    
    figure(k + numel(taps))
    plot(snr_dat(k).freq{1}, squeeze(mean(snr_dat(k).trial(256,:, :), 3)))
    print(gcf, print_format, ['tapers_test_' num2str(taps{k}(1)) '-' num2str(taps{k}(2)) '_snr'])
end
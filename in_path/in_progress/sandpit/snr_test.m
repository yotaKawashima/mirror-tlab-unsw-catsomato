data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr/';
catname = 'C20110808_R03';
area = 'S1';
filename_header = [catname, '_TStim_', area];
print_type = '-dpng';

fname = dir(fullfile(data_dir, [filename_header '*.mat']));

cCond = 16;

load([data_dir, fname(cCond).name])

%lastind = size(data.trial, 2);
lastind = find_closest_ind(data.freq{1}, 300);


%data.trial = mean(data.trial, 3);


snr_out = zeros(size(data.trial, 1), lastind, size(data.trial, 3));

single_f_step = mean(diff(data.freq{1}));
outer = round(3*1/single_f_step);
inner = round(1/single_f_step);

% middle cases
for k = 21:(lastind-21)
    
    snr_out(:, k, :) = data.trial(:, k, :) - mean(data.trial(:, [(k-20):(k-7) (k+7):(k+20)], :), 2);
    
end


% k<21

for k=1:20
    snr_out(:, k, :) = data.trial(:, k, :) - mean(data.trial(:, [1:(k-7) (k+7):(k+20)], :), 2);

end

% k>lastind-21

for k = (lastind-20):(lastind-21)
    
    snr_out(:, k, :) = data.trial(:, k, :) - mean(data.trial(:, [(k-20):(k-7) (k+7):lastind], :), 2);
    
end


%%
snr_plot = mean(snr_out, 3);
snr_plot(isinf(snr_plot))=NaN;
snr_plot = mean(snr_plot(101:end, :), 1, 'omitnan');

figure(1)
plot(data.freq{1}(1:lastind), snr_plot)
shg

%print(1, print_type, 'snr_test_snr')

pwr_plot = mean(data.trial(101:end, 1:lastind, :), 3);
pwr_plot(isinf(pwr_plot))=NaN;
pwr_plot = mean(pwr_plot, 1, 'omitnan');

figure(2)
plot(data.freq{1}(1:lastind), pwr_plot)
shg

%print(2, print_type, 'snr_test_pwr')
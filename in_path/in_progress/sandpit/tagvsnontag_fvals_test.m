
%% Load from file
load('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/compiled_anova_data')
    


%% find the indices of the non-tagged frequencies
[freq_keep, keep_i, hars] = no_resp_indices([23, 200], [], data(1).freq, 'logical');
% data(1).freq == data(2).freq




%% seperate data


for cArea = 1:2
    data(cArea).nontagged.fvals = data(cArea).all.fvals(:, keep_i, :);
    data(cArea).nontagged.freq = freq_keep;
    
    data(cArea).tagged.fvals = data(cArea).all.fvals(:, ~keep_i, :);
    data(cArea).tagged.freq = data(1).freq(~keep_i);
    
end


 save('compiled_anova_data_sepbytag_fvals', 'data')
 
 
 
%% casually combining the data.
data_comb(1).type = 'nontagged';
data_comb(1).fvals = [data(1).nontagged.fvals; data(2).nontagged.fvals];
data_comb(1).freq = data(1).nontagged.freq;
data_comb(1).S1channels = size(data(1).nontagged.fvals, 1);
data_comb(1).S2channels = size(data(2).nontagged.fvals, 1);


data_comb(2).type = 'tagged';
data_comb(2).fvals = [data(1).tagged.fvals; data(2).tagged.fvals];
data_comb(2).freq = data(1).tagged.freq;

%% 

for k = 1:2
    figure(k)
    clf 
    hold on 
    scatter(data_comb(k).freq, permute(mean(data_comb(k).fvals(:, :, 1), 1, 'omitnan'), [2, 3, 1]), 'r.')
    scatter(data_comb(k).freq, permute(mean(data_comb(k).fvals(:, :, 2), 1, 'omitnan'), [2, 3, 1]), 'g.')
    scatter(data_comb(k).freq, permute(mean(data_comb(k).fvals(:, :, 3), 1, 'omitnan'), [2, 3, 1]), 'b.')
    hold off
    shg
end

%%

printtype = '-depsc';
%%
%table_values = zeros(26, 3);
row_headers = [23; 200; 46; 69; 92; 115; 138; 161; 184; 207; ...
    230; 253; 276; 299; 16; 39; 62; 85; 108; 131; 154; 177; 223; 246; 269; 292];

foi_ind = find_closest_ind(data_comb(2).freq, row_headers);

logf = log(data_comb(2).fvals(:, foi_ind, :)+1);

table_values = permute(mean(logf, 1, 'omitnan'), [2, 3, 1]);

fund_inds = [1, 2];
fund_y = table_values(fund_inds, :);
fund_x = row_headers(fund_inds);

hars_inds = [1, 3:14];
hars_y = table_values(hars_inds, :);
hars_x = row_headers(hars_inds);

ints_inds = [2, 15:26];
ints_x = row_headers(ints_inds);
[ints_x, tmpi] = sort(ints_x);
ints_y = table_values(ints_inds(tmpi), :);

figure(4)
plot( hars_x, hars_y, ints_x, ints_y, fund_x, fund_y, 'r*')
legend({'harmonics, cF1', 'harmonics, cF2', 'harmonics, cF1F2', 'intermodulation, cF1', ...
    'intermodulation, cF2', 'intermodulation, cF1F2', 'fundamentals'}, 'Location', 'EastOutside')
title('Trends in f-values')
ylabel('log(f-value + 1)')
xlabel('frequency (Hz)')
print(4, printtype, 'tagged_logfvals_histogram')
shg

figure(5)
histogram(logf)
title('Histogram of f-values')
xlabel('log(f-value + 1)')
ylabel('frequency')
print(5, printtype, 'tagged_logfvals_trends')


%% Figure 4 in paper

[allfreq, sorti] = sort([data_comb(1).freq; data_comb(2).freq]);

space_nontagged = [log(data_comb(1).fvals+1) NaN(2628, length(data_comb(2).freq), 3)];
space_nontagged(:, :, :) = space_nontagged(:, sorti, :);
plot_dat{1} = permute(mean(space_nontagged(1:data_comb(1).S1channels, :, :), 1, 'omitnan'), [2, 3, 1]);
plot_dat{2} = permute(mean(space_nontagged((data_comb(1).S1channels+1):end, :, :), 1, 'omitnan'), [2, 3, 1]);

limits = find_lims(plot_dat);

for k = 1:2
    figure(5+k)
    plot(allfreq, plot_dat{k})
    
    legend({'cF1', 'cF2', 'cF1F2'}, 'Location', 'NorthEast') %, 'Orientation', 'Horizontal')
    xlabel('frequency (Hz)')
    ylabel('f-statistic')
    
    set(gca, 'YLim', limits)
    
    print(5+k, printtype, ['nontagged_logfvals_trends_S' num2str(k)])
end 

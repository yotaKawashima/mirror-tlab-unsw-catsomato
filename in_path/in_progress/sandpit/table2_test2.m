% load data
data_loc = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/';
load([data_loc 'compiled_anova_data_sepbytag.mat'])


% FDR
% use eeglab FDR for the flexibility of para/non-parametric
% downloaded as its own file. 
data_type = {'tagged', 'nontagged'};
for m = 1:numel(data_type)
    % concatenate S1 and S2
    pvals{m} = cat(1, data(1).(data_type{m}).pvals, data(2).(data_type{m}).pvals);
    
    [p_fdr, p_masked] = eeglab_fdr(pvals{m}, 0.05, 'nonparametric');
    
    p_lim(m) = p_fdr;
    p_thresh{m} = p_masked;
    
    p_class1 = double(p_masked);
    p_class1(:, :, 2) = p_class1(:, :, 2)*2;
    p_class1(:, :, 1) = p_class1(:, :, 1)*4;
    p_class1 = sum(p_class1, 3);
    
    p_class{m} = p_class1;
end


% reformat
for m = 1:numel(data_type)
    data_out(m).datatype = data_type{m};
    data_out(m).freq     = data(1).(data_type{m}).freq;
    data_out(m).pvals    = pvals{m};
    data_out(m).plim     = p_lim(m);
    data_out(m).pthresh  = p_thresh{m};
    data_out(m).pclass   = p_class{m};
end

clear pvals p_lim p_thresh p_class

data_out(1).numchanS1 = size(data(1).tagged.pvals, 1);
data_out(1).numchanS2 = size(data(2).tagged.pvals, 1);

data = data_out;

%save('compiled_anova_data_sepbytag_S1S2combined')

%%



table_values = zeros(26, 10);
row_headers = [23; 200; 46; 69; 92; 115; 138; 161; 184; 207; ...
    230; 253; 276; 299; 16; 39; 62; 85; 108; 131; 154; 177; 223; 246; 269; 292];


foi_ind = find_closest_ind(data(1).freq, row_headers);

% standard edges. ensures no bins skipped.
edge_st = -0.5:1:7.5;
% 0 = null
% 1 = interaction only
% 2 = low f dependence only (usually 23). RA-like
% 3 = low f and interaction
% 4 = high f only (usually 200). PC-like
% 5 = high f and interaction
% 6 = low f and high f. linear
% 7 = low f, high f and interaction

for k = 1:numel(foi_ind)
       
    data_tmp = data(1).pclass(:, foi_ind(k));
    
    % counts
    count_tmp = histcounts(data_tmp, edge_st);
    
    table_values(k, [3, 1, 8, 2, 9, 7, 10]) = count_tmp(2:end);
end

numchanstotal = size(data(1).pvals, 1);

table_values(:, 4) = sum(table_values(:, [1, 7, 8, 10]), 2);
table_values(:, 5) = sum(table_values(:, [2, 7, 9, 10]), 2);
table_values(:, 6) = sum(table_values(:, [3, 8, 9, 10]), 2);

table_values = table_values/numchanstotal*100;


%%

fund_inds = [1, 2];
fund_y = table_values(fund_inds, 4:6);
fund_x = row_headers(fund_inds);

hars_inds = [1, 3:14];
hars_y = table_values(hars_inds, 4:6);
hars_x = row_headers(hars_inds);

ints_inds = [2, 15:26];
ints_x = row_headers(ints_inds);
[ints_x, tmpi] = sort(ints_x);
ints_y = table_values(ints_inds(tmpi), 4:6);


plot(fund_x, fund_y, 'r*', hars_x, hars_y, ints_x, ints_y)

legend({'fundamentals', '-', '-', 'harmonics, cF1', 'harmonics, cF2', 'harmonics, cF1F2', ...
    'intermodulation, cF1', 'intermodulation, cF2', 'intermodulation, cF1F2'}, 'Location', 'EastOutside')

shg


%% counts for abstract

% % cF1 at f1, sep by S1/2
% data_tmp = data(1).pclass(:, foi_ind(1));
% S1_data = data_tmp(1:data(1).numchanS1, :);
% S2_data = data_tmp(data(1).numchanS1:end, :);
% 
% count_tmp = histcounts(S1_data, edge_st);
% S1_cF1atf1 = sum(count_tmp([3, 4, 7, 8]));
% count_tmp = histcounts(S2_data, edge_st);
% S2_cF1atf1 = sum(count_tmp([3, 4, 7, 8]));
% 
% 
% % cF2 at f2, sep by S1/2
% data_tmp = data(1).pclass(:, foi_ind(2));
% S1_data = data_tmp(1:data(1).numchanS1, :);
% S2_data = data_tmp(data(1).numchanS1:end, :);
% 
% count_tmp = histcounts(S1_data, edge_st);
% S1_cF2atf2 = sum(count_tmp([5, 6, 7, 8]));
% count_tmp = histcounts(S2_data, edge_st);
% S2_cF2atf2 = sum(count_tmp([5, 6, 7, 8]));




table_values = zeros(26, 10);
row_headers = [23; 200; 46; 69; 92; 115; 138; 161; 184; 207; ...
    230; 253; 276; 299; 16; 39; 62; 85; 108; 131; 154; 177; 223; 246; 269; 292];


foi_ind = find_closest_ind(data(1).freq, row_headers);

% standard edges. ensures no bins skipped.
edge_st = -0.5:1:7.5;
% 0 = null
% 1 = interaction only
% 2 = low f dependence only (usually 23). RA-like
% 3 = low f and interaction
% 4 = high f only (usually 200). PC-like
% 5 = high f and interaction
% 6 = low f and high f. linear
% 7 = low f, high f and interaction

for k = 1:numel(foi_ind)
       
    data_tmp = data(1).pclass(data(1).numchanS1:end, foi_ind(k));
    
    % counts
    count_tmp = histcounts(data_tmp, edge_st);
    
    table_values(k, [3, 1, 8, 2, 9, 7, 10]) = count_tmp(2:end);
end

numchanstotal = size(data(1).pvals, 1);

table_values(:, 4) = sum(table_values(:, [1, 7, 8, 10]), 2);
table_values(:, 5) = sum(table_values(:, [2, 7, 9, 10]), 2);
table_values(:, 6) = sum(table_values(:, [3, 8, 9, 10]), 2);


























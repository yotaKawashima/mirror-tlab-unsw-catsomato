% load data


% count

load('compiled_anova_data_sepbytag_threshold', 'data')

table_values = zeros(26, 36);
row_headers = [23; 200; 46; 69; 92; 115; 138; 161; 184; 207; ...
    230; 253; 276; 299; 16; 39; 62; 85; 108; 131; 154; 177; 223; 246; 269; 292];


foi_ind = find_closest_ind(data(1).tagged.freq, row_headers);

% standard edges. ensures no bins skipped.
edge_st = -0.5:1:7.5;
% 0 = null (we ignore this here)
% 1 = interaction only
% 2 = low f dependence only (usually 23). RA-like
% 3 = low f and interaction
% 4 = high f only (usually 200). PC-like
% 5 = high f and interaction
% 6 = low f and high f. linear
% 7 = low f, high f and interaction

for k = 1:numel(foi_ind)
    
    % -- S1 --
    
    data_tmp = data(1).tagged.thresh(:, foi_ind(k));
    
    % counts
    count_tmp = histcounts(data_tmp, edge_st);
    
    % sort into correct places
    % cF1 only is class = 2, = element 3
    table_values(k, 1) = count_tmp(3);
    
    % cF2 only is 4 = element 5.
    table_values(k, 13) = count_tmp(5);
    
    % cF1F2 only is class 1 = element 2
    table_values(k, 25) = count_tmp(2);
    
    % cF1 in combination is class 2, 3, 6, 7
    table_values(k, 7) = sum(count_tmp([3, 4, 7, 8]));
    
    % cF2 in combination is class 4, 5, 6, 7
    table_values(k, 19) = sum(count_tmp(5:8));
    
    % cF1F2 in combination is class 1,3,5,7
    table_values(k, 31) = sum(count_tmp(2:2:8));
    
    
    % -- S2 --
    
    data_tmp = data(2).tagged.thresh(:, foi_ind(k));
    
    % counts
    count_tmp = histcounts(data_tmp, edge_st);
    
    % sort into correct places
    % cF1 only is class = 2, = element 3
    table_values(k, 3) = count_tmp(3);
    
    % cF2 only is 4 = element 5.
    table_values(k, 15) = count_tmp(5);
    
    % cF1F2 only is class 1 = element 2
    table_values(k, 27) = count_tmp(2);
    
    % cF1 in combination is class 2, 3, 6, 7
    table_values(k, 9) = sum(count_tmp([3, 4, 7, 8]));
    
    % cF2 in combination is class 4, 5, 6, 7
    table_values(k, 21) = sum(count_tmp(5:8));
    
    % cF1F2 in combination is class 1,3,5,7
    table_values(k, 33) = sum(count_tmp(2:2:8));
end
%%
numchanS1 = size(data(1).tagged.pvals, 1);
numchanS2 = size(data(2).tagged.pvals, 1);

table_values(:, [2, 8, 14, 20, 26, 32]) = table_values(:, [1, 7, 13, 19, 25, 31])/numchanS1*100;

table_values(:, [4, 10, 16, 22, 28, 34]) = table_values(:, [3, 9, 15, 21, 27, 33])/numchanS2*100;

numchanstotal = numchanS1 + numchanS2;

table_values(:, 5)  = sum(table_values(:, [1, 3]), 2);
table_values(:, 11) = sum(table_values(:, [7, 9]), 2);
table_values(:, 17) = sum(table_values(:, [13, 15]), 2);
table_values(:, 23) = sum(table_values(:, [19, 21]), 2);
table_values(:, 29) = sum(table_values(:, [25, 27]), 2);
table_values(:, 35) = sum(table_values(:, [31, 33]), 2);

table_values(:, 6:6:36) = table_values(:, 5:6:35)/numchanstotal*100;






































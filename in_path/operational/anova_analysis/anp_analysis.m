function [pvals, tables, stats] = anp_analysis(data_dir, filename_header)

% anp_analysis: analyses data for the anova2 pipeline
%   anp_analysis(data_dir, filename_header) performs anova2 analysis on the
%   data specified by data_dir and filename_header. 
%
%   [pvals, tables, stats] = anp_analysis(...) outputs the analysis to the
%   workspace.

% Written by Rannee Li, Jul 2015
% Change log:
% 10 Jul 2017:  Fixed bug where the name of the file was not updated in
%               metavars.
% 14 Nov 2017:  More general search for filename implemented.
%               Changed load call to be more robust

%% load data
fname = dir(fullfile(data_dir, [filename_header '*_P*anov*adatain.mat']));

if numel(fname) < 1 % more general case
    fname = dir(fullfile(data_dir, [filename_header '*adatain.mat']));
end

load(fullfile(data_dir, fname.name))


%% preallocate variables
nchans = size(data.y, 4);
nfreqs = size(data.y, 3);
pvals = ones(nchans, nfreqs, 3); % holds the p-values
tables = cell(nchans, nfreqs); % holds the data in the table
f = {'source', 'sigmasq', 'colmeans', 'coln', 'rowmeans', 'rown', 'inter', 'pval', 'df'};
stats = struct(f{1}, [], f{2}, [], f{3}, [], f{4}, [], f{5}, [], f{6}, [], f{7}, [], f{8}, [], f{9}, []); % holds variables that can be used later
empt_s = stats;
stats(nchans, nfreqs).df = [];

%% run for each frequency + channel
for m = 1:nfreqs
for k = 1:nchans
    try
        [P,TABLE,STATS] = anova2(data.y(:, :, m, k),data.custom.ntrials,'off');

    catch
        %fprintf('Error for channel %i\n', k)
        P = [NaN, NaN, NaN];
        TABLE = [];
        STATS = empt_s;
    end
    pvals(k, m, :) = P;
    tables{k, m} = TABLE;
    stats(k, m) = STATS;
end
end
metavars = rmfield(data, 'y');
metavars.custom.filename = [fname.name(1:end-11) 'adatout.mat'];

clear data

save(metavars.custom.filename, 'metavars', 'tables', 'pvals', 'stats')

if nargout < 1
    clear pvals tables stats
end

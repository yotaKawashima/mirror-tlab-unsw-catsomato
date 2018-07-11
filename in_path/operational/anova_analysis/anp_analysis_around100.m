function [pvals, tables, stats] = anp_analysis_around100(data_dir, filename_header, FreqSeg)

% anp_analysis: analyses data for the anova2 pipeline
%   anp_analysis(data_dir, filename_header) performs anova2 analysis on the
%   data specified by data_dir and filename_header. 
%
%   [pvals, tables, stats] = anp_analysis(...) outputs the analysis to the
%   workspace.

% Written by Rannee Li, Jul 2015, editted by Phoebe You 19 Feb 2016 for analysing
% interesting response from 50Hz to 150Hz 

%% load data
fname = dir(fullfile(data_dir, [filename_header '*_P*anov*adatain.mat']));
load([data_dir fname.name])

%% Define frequency range
TargFreq = 100;
freq = data.freq{1};
vTargFreq_beforeInterFreq = find(TargFreq-FreqSeg<= freq & TargFreq+FreqSeg>=freq);

% Get rid of 200+n*23 frequecnies
n = -6:1:0;
inter = 200+23.*n;
InterFreq = zeros(size(inter));
% To find the eact match with 200+n*23
for iInter = 1: length(inter)
    InterFreq(iInter) = find_closest (freq, inter(iInter));
end
% Remove intermodutory frequencies
for iInterFreq = 1: length(InterFreq)
    vTargFreq = vTargFreq_beforeInterFreq(freq(vTargFreq_beforeInterFreq)~=InterFreq(iInterFreq));
end

%% preallocate variables
% nchans = size(data.y, 4);
% % nfreqs = size(data.y, 3);
% nfreqs = length(vTargFreq);
% pvals = ones(nchans, nfreqs, 3); % holds the p-values
% tables = cell(nchans, nfreqs); % holds the data in the table
nchans = size(data.y, 4);
ntrials = data.custom.ntrials;
pvals = ones(nchans, ntrials, 3); % holds the p-values
tables = cell(nchans, ntrials); % holds the data in the table

f = {'source', 'sigmasq', 'colmeans', 'coln', 'rowmeans', 'rown', 'inter', 'pval', 'df'};
stats = struct(f{1}, [], f{2}, [], f{3}, [], f{4}, [], f{5}, [], f{6}, [], f{7}, [], f{8}, [], f{9}, []); % holds variables that can be used later
empt_s = stats;
% stats(nchans, nfreqs).df = [];
stats(nchans, ntrials).df = [];

%% run for all channel within frequency range fre(vTargFreq)
% Take mean of power spectrum cross frequency range for each trail each
% channel 
for k = 1:nchans
    for q = 1:ntrials
        try
            [P,TABLE,STATS] = anova2(data.y(:, :, vTargFreq(m), k),data.custom.ntrials,'off');
        catch
        end
    end
end


% for m = 1:nfreqs
% for k = 1:nchans
%     try
%         [P,TABLE,STATS] = anova2(data.y(:, :, vTargFreq(m), k),data.custom.ntrials,'off');
% 
%     catch
%         fprintf('Error for channel %i\n', k)
%         P = [NaN, NaN, NaN];
%         TABLE = [];
%         STATS = empt_s;
%     end
%     pvals(k, m, :) = P;
%     tables{k, m} = TABLE;
%     stats(k, m) = STATS;
% end
% end
metavars = rmfield(data, 'y'); %#ok<NASGU> saved later
clear data

save([fname.name(1:end-11) 'F' num2str(freq(vTargFreq(1))) 'to' num2str(freq(vTargFreq(end))) '_adatout.mat'], 'metavars', 'tables', 'pvals', 'stats')

if nargin < 1
    clear pvals tables stats
end

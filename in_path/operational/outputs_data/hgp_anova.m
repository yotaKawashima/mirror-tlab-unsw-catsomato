function hgp_anova(data, datain_dir)
% load data before using


% preallocate
% the matrix is n23*nTrials x n200 x nChannels elements
nTrials = data.custom.ntrials;
nCol = data.custom.subplotconfig(2);
nRow = data.custom.subplotconfig(1);
nChan = data.custom.nsignals;
y = zeros(nRow*nTrials, nCol, nChan);

% arrange for anova
% skip the first set because it is the baseline and because of the way that
% HGP is calculated from evoked power, there is no basline

for k = 1:numel(data.trial)
    tmpmean = nanmean(data.trial{k}, 2); % mean across frequencies
    
    y(floor(k/nRow)*nTrials+1:(floor(k/nRow)+1)*nTrials, mod(k, nCol)+1, :) = permute(tmpmean, [3, 2, 1]);
    
end

% fix up the output

% variables in the 'main menu'
dat_tmp.conditions  = data.datalabels;
dat_tmp.label       = data.label;
dat_tmp.freq{1}     = data.freq{1};

% data.custom updates
dat_tmp.custom                      = data.custom;
dat_tmp.custom.datatype.isanovadata = true;
dat_tmp.custom                      = rmfield(dat_tmp.custom, {'conditions', 'nsamples'});
dat_tmp.custom.nsignals             = nChan;
dat_tmp.custom.filename             = [data.custom.filename '_adatain'];

% delete comm
data = dat_tmp;
data.y          = y;

% save
save(data.custom.filename, 'data')

afpc_filemover(data.custom.filename(1:13), datain_dir, 'adatain')


% --- and now for the actual anova
nchans = size(data.y, 4);
nfreqs = size(data.y, 3);
pvals = ones(nchans, nfreqs, 3); % holds the p-values
tables = cell(nchans, nfreqs); % holds the data in the table
f = {'source', 'sigmasq', 'colmeans', 'coln', 'rowmeans', 'rown', 'inter', 'pval', 'df'};
stats = struct(f{1}, [], f{2}, [], f{3}, [], f{4}, [], f{5}, [], f{6}, [], f{7}, [], f{8}, [], f{9}, []); % holds variables that can be used later
empt_s = stats;
stats(nchans, nfreqs).df = [];

% run for each frequency + channel
for m = 1:nfreqs
for k = 1:nchans
    try
        [P,TABLE,STATS] = anova2(data.y(:, :, m, k),data.custom.ntrials,'off');

    catch
        fprintf('Error for channel %i\n', k)
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
metavars.custom.filename = [metavars.custom.filename(1:end-8) '_adatout.mat'];

clear data

save(metavars.custom.filename, 'metavars', 'tables', 'pvals', 'stats')

afpc_filemover(metavars.custom.filename(1:13), [datain_dir(1:end-1), '_adatain/'], 'adatout')

if nargout < 1
    clear pvals tables stats
end
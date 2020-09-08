function [data, ttest_data] = f2p_loader(data_dir, save_dir, data_type, max_only, reload, ttest, fmax)
% Figure 2 data creation
%   data_dir: data directory
%   save_dir: where the saved files would be
%   data_type: the folder name/suffix of the required data
%   max_only: true if only the max condition is needed
%   reload: true if the data should be recalculated (usually this does not
%       happen if the data already exists)
%   ttest: true if a ttest should be conducted
%   fmax: the maximum frequency to keep

% get the filename
if max_only
    fname = [data_type '_max.mat'];
else
    fname = [data_type '_all.mat'];
end

% check if files already exist
% r_d :flag for data (true if they already exist.)
% r_t :flag for t-test data (true if they already exist.)
[data, ttest_data, r_d, r_t] = f2pl_check(save_dir, fname, reload);

% if both files are found, then we can just return
if r_d && r_t
    % both data and ttest_data have been found and loaded, so return
    return
end

% find the names of the cats
cat_names = dirsinside(data_dir);

% which files need to be created?
if r_d && ~r_t
    % just need to do the ttest
    
    % check if we need a ttest based on argument of f2p_loader function
    % (input to this whole function.)
    if ttest % 1 when we need ttest_data.
        % load ttest data
        [data, data_full] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
        ttest_data = f2pl_ttest(data, save_dir, fname);
        %ttest_data = f2pl_ttest(data_full, save_dir, fname); %(data_full, true, save_dir, fname);
        return
    else % 0 when we do not need ttest data.
        return
    end
    
elseif ~r_d && r_t
    % just need to load the data
    [data, ~] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
    return
elseif ~r_d && ~r_t
    % need to do everything from fresh
    if ttest 
        [data, data_full] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
        ttest_data = f2pl_ttest(data, save_dir, fname);
        %ttest_data = f2pl_ttest(data_full, save_dir, fname); %(data_full, true, save_dir, fname);
    else
        [data, ~] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
    end
else
    error('Broken logic')    
end


function ttest_data = f2pl_ttest(data, save_dir, fname) %(data, trials, save_dir, fname)
% function to conduct t-test

% number of frequencies
nF = size(data.trial, 2);

% preallocate
hs = zeros(1, nF);
ps = zeros(1, nF);
cis = zeros(2, nF);

% fill the first item
f = 1;
ttestdata = data.trial(:, f); % dim = (channels*sessions*trials) x 1
%{
if trials
    tmp = data.trial(:, 1, :);
    ttestdata = tmp(:);
else
    ttestdata = data.trial(:, 1);
end
%}
ttestdata = ttestdata(~isnan(ttestdata));

[h, p, ci, stats] = ttest(ttestdata, 0, 'Tail', 'right');

stats(nF).tstat = [];

tstats = zeros(1, nF);

hs(f) = h; % 1 if null hypothesis is rejected. 0 if not. (alpha = 5%)
ps(f) = p; % p-value
cis(:, f) = ci; % the lower and upper boundaries of the 100(1-alpha)
% tstat contains tstat, df, and sd.
% tstat - value of the t-statistic
% df - degree of freedom.
% sd - estimated population standard deviation.
tstats(f) = stats.tstat; 


% finish filling everything
for f = 2:nF
    
    % get data
    ttestdata = data.trial(:, f); % dim = (channels*sessions*trials) x 1
    %{
    if trials
        tmp = data.trial(:, f, :);
        ttestdata = tmp(:);
    else
        ttestdata = data.trial(:, f);
    end
    %}
    ttestdata = ttestdata(~isnan(ttestdata));
    
    % conduct t-test
    [h, p, ci, stats] = ttest(ttestdata, 0, 'Tail', 'right');
    
    % put into buckets
    hs(f) = h;
    ps(f) = p;
    cis(:, f) = ci;
    stats(f) = stats;
    tstats(f) = stats.tstat;
end

% arrange
ttest_data.hs = hs;
ttest_data.ps = ps;
ttest_data.cis = cis;
ttest_data.stats = stats;
ttest_data.tstats = tstats;
ttest_data.freq{1} = data.freq{1};

save(fullfile(save_dir, [fname(1:end-4) '_ttest.mat']), 'ttest_data') 




function [data, data_full] = f2pl_load(cat_names, data_dir, data_type, full, max_only, save_dir, fname, fmax)
% data      : mean across trials
%             dim = (channels*sessions) x frequencies
% data_full : data without taking mean across trials
%             dim = (channels*sessions*trials) x frequencies
% preallocate
data_out = [];
data_all = [];
labels = [];
conditions = [];
nsignals = [];
spatialconfigs = [];
areas = [];

for c = 1:numel(cat_names)
    
    % get filenames
    loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));
        
    % check if max condition only
    if max_only 
        nfiles = numel(loadnames);   
        lnames = loadnames;
        % The last file of each area (S1 and S2) is from max condition.
        loadnames = lnames(nfiles/2); % from S1
        loadnames(2) = lnames(nfiles); % from S2
        nfiles = 2; % we have only 2 files. 
    else
        nfiles = numel(loadnames);    
    end
    
    
    for k = 1:nfiles
        % load
        load(fullfile(data_dir, cat_names{c}, data_type, loadnames(k).name))
        
        % work out which channels to keep (the bipolar ones)
        if data.custom.nsignals == 280 || data.custom.nsignals == 176
            lastunipol = prod(data.custom.spatialconfig);
            chan = lastunipol+1:data.custom.nsignals;
        elseif data.custom.nsignals == 180 || data.custom.nsignals == 112
            chan = 1:data.custom.nsignals;
        else
            error('Only unipolar?');
        end
        
        
        if c==1 && k == 1
            [~, f_ind] = find_closest(data.freq{1}, fmax);
        end
        label = data.label(chan);
        labels = [labels; label];
        conditions = [conditions; data.custom.conditions];
        nsignals = [nsignals; data.custom.nsignals];
        spatialconfigs = [spatialconfigs; data.custom.spatialconfig];
        areas = [areas; data.custom.area];
        
        % get it into the avtrcon format
        data_tmp = data.trial(chan, 1:f_ind, :);
        % take mean across trials (repetitions).
        % data_tmp dim = channels x frequencies
        data_tmp = mean(data_tmp, 3); 
        % smush
        % (channels*sessions) x frequencies 
        % Here channels is from S1 and S2.
        data_out = [data_out;data_tmp];
        
        % if the full dataset needed, smush that too
        if full
            if k ~= 1 || k ~= (nfiles/2)+1
                for m = 1:size(data.trial, 3)
                    % (channels*sessions*trials) x frequencies 
                    data_all = [data_all;data.trial(chan, 1:f_ind, m)];
                end
            end
        end
    end
end

% clean data 
data_out(isinf(data_out)) = NaN; % why inf happends?
data.trial = data_out; % dim = (channels*sessions x frequencies)
data.label = labels;
data.custom.filename(2:13) = 'xxxxxxxx_Rxx';
data.freq{1} = data.freq{1}(1:f_ind);
data.label = labels;
data.custom.conditions = conditions;
data.custom.nsignals = nsignals;
data.custom.spatialconfigs = spatialconfigs;
data.custom.areas = areas;

save(fullfile(save_dir, fname), 'data') 


% clean data for full 
if full
    data_all(isinf(data_all)) = NaN;
    data_full = data;
    data_full.trial = data_all;
    data_full.custom.filename(2:13) = 'xxxxxxxx_Rxx';
    data_full.freq{1} = data.freq{1}(1:f_ind);
    data_full.label = labels;
    data_full.custom.conditions = conditions;
    data_full.custom.nsignals = nsignals;
    data_full.custom.spatialconfigs = spatialconfigs;
    data_full.custom.areas = areas;
    
else
    data_full = [];
end

function [data, ttest_data, r_d, r_t] = f2pl_check(save_dir, fname, reload)

% set defaults
r_d = false;
r_t = false;
data = [];
ttest_data = [];


% check if data file already exists
name = dir(fullfile(save_dir, fname));
if ~isempty(name) && ~reload
    load(fullfile(save_dir, name.name))
    r_d = true;
end

% check if ttest file already exists
name = dir(fullfile(save_dir, [fname(1:end-4) '_ttest.mat']));
if ~isempty(name) && ~reload
    load(fullfile(save_dir, name.name))
    r_t = true;
end













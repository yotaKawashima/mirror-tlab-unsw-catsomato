function [data, ttest_data] = f2p_loader(data_dir, save_dir, data_type, max_only, reload, ttest, fmax)

if max_only
    fname = [data_type '_max.mat'];
else
    fname = [data_type '_all.mat'];
end

% check if files already exist
[data, ttest_data, r_d, r_t] = f2pl_check(save_dir, fname, reload);



if r_d && r_t
    % both data and ttest_data have been found and loaded, so return
    return
end

% find the names of the cats
cat_names = dirsinside(data_dir);

if r_d && ~r_t
    % just need to do the ttest
    
    % check if we need a ttest
    if ttest
        [~, data_full] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
        ttest_data = f2pl_ttest(data_full, true, save_dir, fname);
        return
    else
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
        ttest_data = f2pl_ttest(data_full, true, save_dir, fname);
    else
        [data, ~] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax);
    end
else
    error('Broken logic')    
end


function ttest_data = f2pl_ttest(data, trials, save_dir, fname)
if nargin < 2
    trials = true;
end

nF = size(data.trial, 2);

hs = zeros(1, nF);
ps = zeros(1, nF);
cis = zeros(2, nF);

f = 1;

if trials
    tmp = data.trial(:, 1, :);
    ttestdata = tmp(:);
else
    ttestdata = data.trial(:, 1);
end

[h, p, ci, stats] = ttest(ttestdata, 0, 'Tail', 'right');

stats(nF).tstat = [];

tstats = zeros(1, nF);

hs(f) = h;
ps(f) = p;
cis(:, f) = ci;
tstats(f) = stats.tstat;

for f = 2:nF
    if trials
        tmp = data.trial(:, f, :);
        ttestdata = tmp(:);
    else
        ttestdata = data.trial(:, f);
    end
    
    ttestdata = ttestdata(~isnan(ttestdata));

    [h, p, ci, stats] = ttest(ttestdata, 0, 'Tail', 'right');

    hs(f) = h;
    ps(f) = p;
    cis(:, f) = ci;
    stats(f) = stats;
    tstats(f) = stats.tstat;
end

ttest_data.hs = hs;
ttest_data.ps = ps;
ttest_data.cis = cis;
ttest_data.stats = stats;
ttest_data.tstats = tstats;
ttest_data.freq{1} = data.freq{1};

save(fullfile(save_dir, [fname(1:end-4) '_ttest.mat']), 'ttest_data') 




function [data, data_full] = f2pl_load(cat_names, data_dir, data_type, full, max_only, save_dir, fname, fmax)
data_out = [];
data_all = [];

for c = 1:numel(cat_names)
    % get filenames
    loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));
    
    % check if max condition only
    nfiles = numel(loadnames);
    if max_only
        lnames = loadnames;
        loadnames = lnames(nfiles/2);
        loadnames(2) = lnames(nfiles);
        nfiles = 2;
    end
    
    
    for k = 1:nfiles
        % load
        load(fullfile(data_dir, cat_names{c}, data_type, loadnames(k).name))
        
        if c==1 && k == 1
            [~, f_ind] = find_closest(data.freq{1}, fmax);
        end
        
        % get it into the avtrcon format
        data_tmp = data.trial(:, 1:f_ind, :);
        data_tmp = mean(data_tmp, 3);
        
        % smush
        data_out = [data_out;data_tmp];
        if full
            if k ~= 1 || k ~= (nfiles/2)+1
                for m = size(data.trial, 3)
                    data_all = [data_all;data.trial(:, 1:f_ind, m)];
                end
            end
        end
    end
end

data_out(isinf(data_out)) = NaN;
data.trial = nanmean(data_out, 1);
data.custom.filename(2:13) = 'xxxxxxxx_Rxx';
data.freq{1} = data.freq{1}(1:f_ind);

save(fullfile(save_dir, fname), 'data') 

if full
    data_all(isinf(data_all)) = NaN;
    data_full = data;
    data_full.trial = data_all;
    data_full.custom.filename(2:13) = 'xxxxxxxx_Rxx';
    data.freq{1} = data.freq{1}(1:f_ind);
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













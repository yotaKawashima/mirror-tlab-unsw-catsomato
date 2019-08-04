function [data, ttest_data] = f2p_loader(data_dir, save_dir, data_type, max_only, reload, ttest, fmax, cat_name, channel)
% Figure 2 data creation
%   data_dir: data directory
%   save_dir: where the saved files would be
%   data_type: the folder name/suffix of the required data
%   max_only: true if only the max condition is needed
%   reload: true if the data should be recalculated (usually this does not
%       happen if the data already exists)
%   ttest: true if a ttest should be conducted
%   fmax: the maximum frequency to keep
%   cat_name: name of the cat files. set to false to use all cats
%   channel: the channel of data to load. set to <= 0 to use all channels

% get the filename
if max_only
    fname = [data_type '_max.mat'];
else
    fname = [data_type '_all.mat'];
end

% check if files already exist
[data, ttest_data, r_d, r_t] = f2pl_check(save_dir, fname, reload);

% if both files are found, then we can just return
if r_d && r_t
    % both data and ttest_data have been found and loaded, so return
    return
end

% find the names of the cats
if not cat_name
    cat_names = dirsinside(data_dir);
else 
    cat_names = {cat_name};

% which files need to be created?
if r_d && ~r_t
    % just need to do the ttest
    
    % check if we need a ttest
    if ttest
        [~, data_full] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax, channel);
        ttest_data = f2pl_ttest(data_full, true, save_dir, fname);
        return
    else
        return
    end
    
elseif ~r_d && r_t
    % just need to load the data
    [data, ~] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax, channel);
    return
elseif ~r_d && ~r_t
    % need to do everything from fresh
    if ttest
        [data, data_full] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax, channel);
        ttest_data = f2pl_ttest(data_full, true, save_dir, fname);
    else
        [data, ~] = f2pl_load(cat_names, data_dir, data_type, ttest, max_only, save_dir, fname, fmax, channel);
    end
else
    error('Broken logic')    
end


function ttest_data = f2pl_ttest(data, trials, save_dir, fname)
% function to conduct t-test

% number of frequencies
nF = size(data.trial, 2);

% preallocate
hs = zeros(1, nF);
ps = zeros(1, nF);
cis = zeros(2, nF);

% fill the first item
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

% finish filling everything
for f = 2:nF
    
    % get data
    if trials
        tmp = data.trial(:, f, :);
        ttestdata = tmp(:);
    else
        ttestdata = data.trial(:, f);
    end
    
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




function [data, data_full] = f2pl_load(cat_names, data_dir, data_type, full, max_only, save_dir, fname, fmax, channel)

% preallocate
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

            % finally deal with the channel
            if channel > 0
                c_ind = channel;% TODO
            else
                c_ind = 1:len(data.trial, 1)
            end

        end
        
        % get it into the avtrcon format
        data_tmp = data.trial(c_ind, 1:f_ind, :);
        data_tmp = mean(data_tmp, 3);
        
        % smush
        data_out = [data_out;data_tmp];
        
        % if the full dataset needed, smush that too
        if full
            if k ~= 1 || k ~= (nfiles/2)+1
                for m = 1:size(data.trial, 3)
                    data_all = [data_all;data.trial(c_ind, 1:f_ind, m)];
                end
            end
        end
    end
end

% clean data
data_out(isinf(data_out)) = NaN;
data.trial = nanmean(data_out, 1);
data.custom.filename(2:13) = 'xxxxxxxx_Rxx';
data.freq{1} = data.freq{1}(1:f_ind);
data.custom.channels_out = c_ind;

save(fullfile(save_dir, fname), 'data') 


% clean data for full 
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













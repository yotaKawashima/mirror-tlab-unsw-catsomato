function data = collectallcats(data_dir, data_type, varargin)

% collectallcats: collates data scross cats produced by collectallchannels
%
% collectallcats(data_dir, data_type) collects the data as specified by
% the path [data_dir, data_type]. The output is saved to file.
%
% collectallcats(data_dir, data_type, 'name', arg ...) takes in the
% name-argument pairs as specified below.
%   fromfile: str of the file to load, if it already exists. Suppresses
%       data collation. Default [].
%   ifttest: logical. If true, performs a t-test. Forced to be true if
%       plotttest is set to a non-false value. Default false.
%   plotttest: logical or numeric. If true, plots the t-stat on a new
%       figure. If numeric, plots on that figure. Default false.
%   plotmean: logical or numeric. Same as plotttest, but for mean across
%   	channels. Default false.
%   dontsave: logical or numeric. overwrites default save behaviour. When
%       true, does not save the data. 
%   domain
%
% data = collectallcats(...) does not save the data to a file and instead
% outputs it to the workspace in the struct data.


% % data directory
% dir_head = '/media/rannee/UNSW_Cat_Somatos/';
% data_dir = [dir_head '/data/collated_data/'];
% data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_avtrcon/';


%% input parsing
p = inputParser;
addRequired(p, 'data_dir', @ischar);
addRequired(p, 'data_type', @ischar);

charoremp = @(x) ischar(x) || isempty(x);
addParameter(p, 'fromfile', [], charoremp); % load data from file

numorlog = @(x) isnumeric(x) || islogical(x);
addParameter(p, 'ifttest', false, numorlog);
addParameter(p, 'plotttest', false, numorlog);
addParameter(p, 'plotmean', false, numorlog);
addParameter(p, 'trials', false, numorlog);
addParameter(p, 'dontsave', false, numorlog);
addParameter(p, 'domain', [], @isnumeric);

parse(p, data_dir, data_type, varargin{:})

fromfile    = p.Results.fromfile;
ifttest     = p.Results.ifttest;
plotmean    = p.Results.plotmean;
plotttest   = p.Results.plotttest;
trials      = p.Results.trials;
dontsave    = p.Results.dontsave;
domain      = p.Results.domain;

%% load, collate and save
if numel(fromfile)>0
    % load existing file
    load(fullfile(data_dir, fromfile))
    
elseif strcmp(data_type(end-6:end), 'avtrcon')
    % make file
    
    % find cats
    loadname = dir(fullfile(data_dir, data_type, '*.mat'));
    
    data_out = [];
    
    for c = 1:numel(loadname)
        % load
        load(fullfile(data_dir, data_type, loadname(c).name))
        
        % smush
        data_out = [data_out;data.trial];
        
    end
    
    data_out(isinf(data_out)) = NaN;
    
    data.trial = data_out;
    data.custom.filename(2:13) = 'xxxxxxxx_Rxx';
    
    if nargout < 1 && ~dontsave
        % save this data
        save(data.custom.filename, 'data')
    end
else % generic condition
    
    data_out = [];
    
    cat_names = dirsinside(data_dir);
    
    for c = 1:numel(cat_names)
        
        loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));
        
        for k = 1:numel(loadnames)
            % load
            load(fullfile(data_dir, cat_names{c}, data_type, loadnames(k).name))

            % get it into the avtrcon format
            data_tmp = data.trial;
            data_tmp = mean(data.trial, 3);

            % smush
            data_out = [data_out;data_tmp];
        end
    end
        
    data_out(isinf(data_out)) = NaN;
    data.trial = data_out;
    data.custom.filename(2:13) = 'xxxxxxxx_Rxx';
    
    if nargout < 1 && ~dontsave
        % save this data
        save(data.custom.filename, 'data')
    end
end

%% plot mean

if plotmean
    if islogical(plotmean)
        figure()
    else
        figure(plotmean)
    end
    
    if isempty(domain)
        plot(data.freq{1}, nanmean(data.trial, 1))
    else
        [~, a] = find_closest(data.freq{1}, domain(1));
        [~, b] = find_closest(data.freq{1}, domain(end));
        plot(data.freq{1}(a:b), nanmean(data.trial(:,a:b,:), 1))
    end
    
    xlabel('frequency')
    ylabel('power')
    
end

%% do the t-test, if required.
if ifttest || (~ifttest&&plotttest)
        
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
        
        [h, p, ci, stats] = ttest(ttestdata, 0, 'Tail', 'right');
        
        hs(f) = h;
        ps(f) = p;
        cis(:, f) = ci;
        stats(f) = stats;
        tstats(f) = stats.tstat;
    end
    
    if nargout < 1
        save([data.custom.filename '_ttested'], 'hs', 'ps', 'cis', 'stats', 'tstats')
    end
end


if plotttest
    if islogical(plotttest)
        figure()
    else
        figure(plotttest)
    end
    
    plot(data.freq{1}, tstats)
    
    xlabel('frequency')
    ylabel('t-statistic')
end

%% deal with output
if nargout < 1
    clear data
else
    data.ttest.hs = hs;
    data.ttest.ps = ps;
    data.ttest.cis = cis;
    data.ttest.stats = stats;
    data.ttest.tstats = tstats;
end


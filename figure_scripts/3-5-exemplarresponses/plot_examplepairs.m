% Plots example response

%% Data selection
cat_name = 'C20110808_R03';

% area, channel, frequency
plotpair{1} = [2, 49, 200];
plotpair{2} = [1, 51, 23];
plotpair{3} = [1, 33, 200];

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%
%img_fmt = '-dpng';
img_fmt = '-dsvg';

%% Set script specific variables
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';
p_data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr';

%% Find files

% find list of names
data_dir = fullfile(data_path, 'included_datasets');
loadnames = dir(fullfile(data_dir, cat_name, data_type, '*.mat'));

% check if the data files exist! if not, make them.
if isempty(loadnames)
    call_snrtosurrounds(data_path, data_type(1:end-8), cat_name)
    
    loadnames = dir(fullfile(data_path, cat_name, data_type, '*.mat'));
end

% find some numbers
nCond = numel(loadnames)/2; % number of conditions per area 
allps = cell(1, 2); % one for each area
allhs = cell(1, 2);
meansnrs = cell(1, 2);
snrline = cell(1, 2);

%% Parse input
FOIs = cell(2,1); % 2 comes from number of area (S1 and S2).
chans = cell(2,1); % Store information for each area.
for p = 1:numel(plotpair)
    % extract the frequencies of interest
    FOIs{plotpair{p}(1)} = [FOIs{plotpair{p}(1)}, plotpair{p}(3)];
    chans{plotpair{p}(1)} =[chans{plotpair{p}(1)}, plotpair{p}(2)];
end

%% Get data

for a = 1:2
    FOI = FOIs{a};
    nFOI = numel(FOI);
    
    for c = 1:nCond
        
        % load
        fprintf('Loading area %i/%i, loading data %2i/%i\n', a, 2, c, nCond)
        load(fullfile(data_dir, cat_name, data_type, loadnames((a-1)*nCond + c).name))
        
        % if it's the first, check whether data included unipolar and
        % preallocate frequencies of interest (per area).
        if c==1
            % check if unipolar id included in srn data. If included,
            % unipolar should be removed. 
            if strcmp(data.label{1}(1:3), 'raw')
                c1_ind = 1 + prod(data.custom.spatialconfig); % the first bipolar channel
                nChan = data.custom.nsignals - c1_ind + 1; % # of bipolar channels
            else
                nChan = data.custom.nsignals;
                c1_ind = 1;
            end
            
            % grab data
            allps{a} = zeros(nChan, nFOI, nCond); % bipolar ch x target x conditions
            allhs{a} = zeros(nChan, nFOI, nCond);
            meansnrs{a} = zeros(nChan, nFOI, nCond);
                        
            iFOI = zeros(1, nFOI); % id of freq
            vFOI = zeros(1, nFOI); % value of freq
            for k = 1:nFOI
                [vFOI(k), iFOI(k)] = find_closest(data.freq{1}, FOI(k));
            end
            
        else % c>1
            % check that the value of frequency is right
            if ~isequal(data.freq{1}(iFOI), vFOI')
                warning('Frequency indicies inconsistent. Recalculating. (a=%i, c=%i)', a, c)
                
                for k = 1:nFOI
                    [vFOI(k), iFOI(k)] = find_closest(data.freq{1}, FOI(k));
                end
            end
        end % c=1:nCond
        
        % extract snrs
        snrs = data.trial(c1_ind:end, iFOI, :);
                
        % do the t-test per channel and per frequency 
        for k = 1:nChan % channel
            for f = 1:nFOI % frequency
                [allhs{a}(k,f,c) , allps{a}(k,f,c)]=ttest(snrs(k, f, :), 0, 'Tail', 'right');
                
            end
        end
        
        % find mean across trials per channel per frequency per area per
        % condition.
        meansnrs{a}(:, :, c) = mean(snrs, 3);
        
        
        % grab the raw data for plotting too if data_type and p_data_type
        % are the same.
        if strcmp(p_data_type, data_type)
            if c==1
                snrline{a} = zeros(nChan, size(data.trial, 2), nCond);
            end
            snrline{a}(:, :, c) = mean(data.trial(c1_ind:end, :, :), 3);
        end
    end
    
    labels{a} = data.label(c1_ind:end);
    spatialconfigs{a} = data.custom.spatialconfig; % channel config
    condfigs{a} = data.custom.subplotconfig; % condition config
end

%% Load plot data id needed
if ~strcmp(p_data_type, data_type)
    % if the data type for the plot is not the same as the spatial
    % distribution
    
    p_data_dir = fullfile(data_path, 'included_datasets');
    p_loadnames = dir(fullfile(p_data_dir, cat_name, p_data_type, '*.mat'));
    
    for a = 1:2
        FOI = FOIs{a};
        nFOI = numel(FOI);
        
        for c = 1:nCond
            
            % load
            fprintf('Loading area %i/%i, loading data %2i/%i\n', a, 2, c, nCond)
            load(fullfile(p_data_dir, cat_name, p_data_type, p_loadnames((a-1)*nCond + c).name))
            
            if c==1
                % check if unipolar id included in srn data. If included,
                % unipolar should be removed. 
                if strcmp(data.label{1}(1:3), 'raw')
                    c1_ind = 1 + prod(data.custom.spatialconfig);
                    nChan = data.custom.nsignals - c1_ind + 1;
                else
                    nChan = data.custom.nsignals;
                    c1_ind = 1;
                end
                snrline{a} = zeros(nChan, size(data.trial, 2), nCond);
            end
            
            snrline{a}(:, :, c) = mean(data.trial(c1_ind:end, :, :), 3);
        end
    end
    
end

%% Threshold SNRs
q = 0.05;

for a = 1:2
    meansnrs{a}(allhs{a}==0)=NaN;
    
    FOI = FOIs{a};
    nFOI = numel(FOI);
    
    for f = 1:nFOI
        % do FDR
        % pID is a threshold for p-val. (scalar not vector)
        [pID, ~] = eeglab_fdr(allps{a}(:, f, :), q, 'parametric');
        tmp = meansnrs{a}(:, f, :); % channel x 1 (f) x conditions.
        tmp(allps{a}(:, f, :)>=pID)=NaN; % Set NaN for NOT significant ones.
        
        meansnrs{a}(:, f, :) = tmp;
        
        fprintf('S%i f=%3i: %f\n', a, FOI(f), pID)
    end % f = 1:nFOI frequency (23 or 200)
    
end % a = 1:2 area

%% Plot spatial map
fnum=1; % first figure number


nAreas = numel(meansnrs);
[nChans, nFOI, nConds] = size(meansnrs{1});

clim = find_lims(meansnrs);
cmap = [0.8275*[1 1 1]; hot(ceil(diff(clim))*10)];

for a = 1:nAreas
    
    FOI = FOIs{a};
    nFOI = numel(FOI);
    
    % relabel channels
    relabels = draw_biprref_chlabfunction(labels{a});
    condfig = condfigs{a};
    
    for f = 1:nFOI
        fig = (a-1)*nAreas + (f-1) + fnum; % has bugs if too high but works for now
        % make a figure
        figure(fig)
        clf
        set(gcf, 'Name', ['f' num2str(FOI(f)) ' S' num2str(a)])
        
        % find channel of interest location
        channo = chans{a}(f);
        [x1, y1] = patch_helper(channo, relabels);
        
        for co = 1:nConds
            % make subplot
            subtightplot(condfig(1), condfig(2), co)
            % goes across rows
            
            values = meansnrs{a}(:,f,co);
            [p, FV] = draw_biprref(values, relabels, spatialconfigs{a}, clim, true);
            colormap(cmap)
            
            % add patch to foi
            p1 = patch(x1, y1, 12);
            set(p1, 'EdgeColor', 'g', 'LineWidth', 2, 'FaceColor', 'none')
        end
        
        % call function to clean up and label the plot
        subtightplotcleaner(fig, condfig, 'cleanticks', false, ...
            'catnames', {loadnames((a-1)*nConds+1:a*nConds).name}, ...
            'topinds', 27:34, 'sideinds', 18:25, 'box', false)
        
        % save figure
        print(gcf, img_fmt, ['Fig4-5_f' num2str(FOI(f)) '_S' num2str(a)])
    end
end


% Save a color bar image.
figure((a-1)*nAreas + (f-1) + fnum+1)
colorbar
colormap(cmap)
caxis(clim)
print(gcf, img_fmt, 'Fig4-5_colorbar')

%% Plot single channels
fnum = (a-1)*nAreas + (f-1) + fnum+2;
% plotpair{1} = [2, 49, 200];

for a = 1:nAreas
    
    FOI = FOIs{a};
    nFOI = numel(FOI);
    
    % relabel channels
    relabels = draw_biprref_chlabfunction(labels{a});
    condfig = condfigs{a};
    
    for f = 1:nFOI
        fig = (a-1)*nAreas + (f-1) + fnum; % has bugs if too high but works for now
        % make a figure
        figure(fig)
        clf
        set(gcf, 'Name', ['f' num2str(FOI(f)) ' S' num2str(a)]) 
        
        
        % frequency constants
        [~, f_i] = find_closest(data.freq{1}, FOI(f));
        iss = mean(diff(data.freq{1}));
        hbw_inds = floor(3/iss);
        freqs = data.freq{1}(f_i-hbw_inds:f_i+hbw_inds);
        
        vals = snrline{a}(chans{a}(f),f_i-hbw_inds:f_i+hbw_inds,:);
        lims = find_lims(vals);
        lims(1) = lims(1) - 1;
        lims(2) = lims(2) + 1;
        
        for co = 1:nConds
            % make subplot
            subtightplot(condfig(1), condfig(2), co)
            % goes across rows
            
            values = vals(:, :, co);
            plot(freqs, values)
            
            %[p, FV] = draw_biprref(values, relabels, spatialconfigs{a}, clim, true);
            %colormap(cmap)
        end
        
        % call function to clean up and label the plot
        subtightplotcleaner(fig, condfig, 'cleanticks', true, ...
            'catnames', {loadnames((a-1)*nConds+1:a*nConds).name}, ...
            'topinds', 27:34, 'sideinds', 18:25, 'box', true, ...
            'xaxisscale', [freqs(1), freqs(end)], 'yaxisscale', lims)
        
        % save figure
        print(gcf, img_fmt, ['Fig4-5_f' num2str(FOI(f)) '_S' num2str(a) '_ch' num2str(chans{a}(f))])
    end
end

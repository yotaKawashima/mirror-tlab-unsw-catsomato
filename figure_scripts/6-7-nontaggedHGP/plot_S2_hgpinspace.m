% Plots spacial distribution of channels with mean high gamma power

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Script specific variables
cat_name = 'C20110808_R03';
ch = 102;
data_dir = [data_path 'included_datasets/' cat_name '/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/'];

filename = [cat_name '_S2_FxxxAxxx_FxxxAxxx_epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp.mat'];

%img_fmt = '-depsc';
img_fmt = '-dsvg';

%% load data 
load([data_dir, filename])

%% preproc
% work out which channels to keep (the bipolar ones)
if strcmp(data.label{1}(1:3), 'raw')
    c1_ind = 1 + prod(data.custom.spatialconfig);
    nChan = data.custom.nsignals - c1_ind + 1;
else
    nChan = data.custom.nsignals;
    c1_ind = 1;
end

% frequency is already selected

% grab extra variables
labels = data.label(c1_ind:end);
spatialconfigs = data.custom.spatialconfig;
condfigs = data.custom.subplotconfig;

%% loop through conditions
for t = 1:numel(data.trial)
    %% select the right channels
    hgp_power = mean(data.trial{t}(c1_ind:end, :, :), 2);
    
    %% perform t-test
    for k = 1:nChan
        [allhs(k,1,t) , allps(k,1,t)]=ttest(hgp_power(k, 1, :), 0, 'Tail', 'right');
    end
    
    %% find means
    mean_hgp(:, :, t) = mean(hgp_power, 3);
end

%% threshold
q = 0.05;

mean_hgp(allhs==0)=NaN;

% do FDR
[pID, ~] = eeglab_fdr(allps(:, 1, :), q, 'parametric');
tmp = mean_hgp(:, 1, :);
tmp(allps(:, 1, :)>=pID)=NaN;

mean_hgp(:, 1, :) = tmp;

%% plot

clim = find_lims(mean_hgp);
clim(2) = clim(2);

cmap = [0.8275*[1 1 1]; hot(ceil(diff(clim))*10)];

nAreas = numel(mean_hgp);
[nChans, ~, nConds] = size(mean_hgp);

% relabel channels
relabels = draw_biprref_chlabfunction(labels);

% make a figure
figure(1)
clf
set(gcf, 'Name', 'mean HGP S2')

[x1, y1] = patch_helper(ch, relabels);

for co = 1:nConds
    % make subplot
    subtightplot(condfigs(1), condfigs(2), co+1)
    % goes across rows
    
    draw_biprref(mean_hgp(:,1,co), relabels, spatialconfigs, clim, true)
    colormap(cmap)
    
    % add patch to foi
    p1 = patch(x1, y1, 12);
    set(p1, 'EdgeColor', 'g', 'LineWidth', 2, 'FaceColor', 'none')
    
end

% call function to clean up and label the plot
subtightplotcleaner(1, condfigs, 'cleanticks', false, ...
    'catnames', [{'F023A000_F200A000'}, data.datalabels], ...
    'topinds', 10:17, 'sideinds', 1:8, 'box', false)

dt = datestr(now, 'yyyymmdd');

% save figure
print(gcf, img_fmt, ['Fig6d_S2_' dt])


% Separately save color bar for this plot 
figure(2)
colorbar
colormap(cmap)
caxis(clim)
print(gcf, img_fmt, ['Fig6d_colorbar_' dt])

% Plots the evoked power in each condition across 50-150 Hz

%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% File specific variables
cat_name = 'C20110808_R03';
ch = 102;

data_dir = [data_path 'included_datasets/' cat_name '/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/'];

filename = [cat_name '_S2_FxxxAxxx_FxxxAxxx_epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp.mat'];

img_fmt = '-depsc';
%% Load and pre-process

load([data_dir, filename])

% Work out if monopolar channels are still present in the data
if strcmp(data.label{1}(1:3), 'raw')
    % monopolar channels present. remove them. 
    c1 = prod(data.custom.spatialconfig)+1;
    data.label = data.label(c1:end);
    for k = 1:numel(data.trial)
        data.trial{k} = data.trial{k}(c1:end, :, :);
    end
end

% Get some required variables
rows = data.custom.subplotconfig(1);
cols = data.custom.subplotconfig(2);
nCond = rows*cols;
plotdata = data.trial;

%% Check if the channel is suitable
% perform a t-test
[h, p] = ttest(mean(plotdata{nCond-1}(ch, :, :), 2), 0, 'Tail', 'right');

% If channel is not significant, don't bother plotting.
if h == 0 || p>0.05
    fprintf('Skipping channel %i (%s)\n', ch, data.label{ch})
    return
else
    fprintf('t-test passed. Plotting channel %i (%s)\n', ch, data.label{ch})
end


%% Plot
figure(ch); clf

yls = zeros(nCond-1, 2);

for p = 1:nCond-1
    subtightplot(rows, cols, p+1)
    plot(data.freq{1}, mean(plotdata{p}(ch,:,:),3))
    yls(p, :) = find_lims(mean(plotdata{p}(ch,:,:),3));
end

ylims = find_lims(yls);
ylims(1) = ylims(1)-1;
ylims(2) = ylims(2)+1;

catnames = cell(numel(data.trial)+1, 1);
catnames{1} = 'F023A000_F200A000';
catnames(2:numel(data.trial)+1) = data.datalabels;

subtightplotcleaner(ch, [rows, cols], 'catnames', catnames, ...
    'xaxisscale', 'tight', 'yaxisscale', ylims, ...
    'sideinds', 1:8, 'topinds', 10:17)


str = sprintf('Fig6c_%s_S2_Ch%03i_highgammapower', cat_name, ch);
print(gcf, img_fmt, str);
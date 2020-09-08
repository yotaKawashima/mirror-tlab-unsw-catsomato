% Plots the evoked power in each condition across 50-150 Hz
% The condition (F1 amplitude ,F2 amplitude) = (0, 0) is included.
%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%% File specific variables
cat_name = 'C20110808_R03';
ch = 102;

data_dir = [data_path 'included_datasets/' cat_name '/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/'];

filename = [cat_name '_S2_FxxxAxxx_FxxxAxxx_epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp.mat'];

img_fmt = '-depsc';
%img_fmt = '-dpng';

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

%% Plot
%for ch = 1:size(plotdata{1}, 1)
f = figure('visible', 'off'); clf

yls = zeros(nCond-1, 2);

for condid = 1:nCond
    % Subplot.
    hAx(condid) = subplot(4, 4, condid);
    % Plot VELogP in 50Hz<f<150Hz. Frequencies of interest (i.e. tagged
    % frequencies) are not included.
    % Take mean and std across trials.
    mean_ = mean(plotdata{condid}(ch, :, :), 3); 
    std_ = std(plotdata{condid}(ch, :, :), 0, 3);
    upper = mean_ + std_;
    lower = mean_ - std_;
    hold on
    % Show std with gray shading
    h = fill([data.freq{1}, fliplr(data.freq{1})], ...
        [upper, fliplr(lower)], 'r');
    set(h, 'FaceColor', [192,192,192]/256);
    set(h, 'EdgeColor', [192,192,192]/256);

    % Show mean with a black line
    plot(data.freq{1}, mean_, 'color', 'k', 'LineWidth', 1);
    hold off
    yls(condid, :) = [min(lower), max(upper)];

end

% Get min and max across stimulus conditions 
ylims = find_lims(yls);
%ylims(1) = ylims(1)-1;
ylims(1) = -10;
ylims(2) = ylims(2)+1;

% Set x and y axis
xlim(hAx, [50 150]);
ylim(hAx, ylims);

% save figure
catnames = cell(numel(data.trial)+1, 1);
catnames{1} = 'F023A000_F200A000';
catnames(2:numel(data.trial)+1) = data.datalabels;

str = sprintf('Fig9d_%s_S2_Ch%03i_highgammapower', cat_name, ch);
print(gcf, img_fmt, str);
close gcf;
%end % for ch = 1:size(plotdata{1}, 1)

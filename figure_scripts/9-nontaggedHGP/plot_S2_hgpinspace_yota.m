% Plots spacial distribution of channels with mean high gamma power
% The condition (F1 amplitude ,F2 amplitude) = (0, 0) is included.
%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Script specific variables
cat_name = 'C20110808_R03';
ch = 102; % The channels chosen in S2_hgp_powerbycond_yota.m
data_dir = [data_path 'included_datasets/' cat_name '/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/'];

filename = [cat_name '_S2_FxxxAxxx_FxxxAxxx_epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp.mat'];

img_fmt = '-depsc';
%img_fmt = '-dpng';

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
for condition = 1:numel(data.trial)
    % select the right channels
    % mean across freq as we difine HGP
    hgp_power = mean(data.trial{condition}(c1_ind:end, :, :), 2); 
    
    % find means across trials within the same condition.
    % bipolar channels x conditions
    mean_hgp(:, condition) = squeeze(mean(hgp_power, 3));
end

%% plot spatial map
% Get parameters for plot
spatialconfig = data.custom.spatialconfig; % Bipolar channel info
condconfig = data.custom.subplotconfig; % Stimulus condition
[nChans, nConds] = size(mean_hgp);

% relabel channels
relabels = draw_biprref_chlabfunction(labels);

% Get exemplar bipolar channel location
[x1, y1] = patch_helper(ch, relabels);

% make a figure
figure(1);clf
set(gcf, 'Name', 'mean HGP S2')

% Preallocate subplots.
hAx = gobjects(nConds,1);

% Remove inf and set color lim if exists
mean_hgp_ = mean_hgp;
if isinf(mean_hgp_)
    mean_hgp_(isinf(mean_hgp_)) = NaN;
    warning('inf is detected. Some bipolar channles may be invalid.');
end %isinf(mean_hgp_)
clim = [0 max(nanmax(mean_hgp_))];
colormap(flipud(hot));

for condid = 1:nConds
    % Set subplot
    hAx(condid) = subtightplot(condconfig(1), condconfig(2), condid);

    % Goes across rows
    values = mean_hgp_(:, condid);
    draw_biprref_vy(values, relabels, spatialconfig, clim);
    hold on
    % add patch to an exemplar bipolar channel
    p1 = patch(x1, y1, 12);
    set(p1, 'EdgeColor', [0,0,125]/256, 'LineWidth', 2, 'FaceColor', 'none')
    hold off
end

% Make axsis invisible
set(hAx, 'XTicklabel',[])
set(hAx, 'YTicklabel',[])

% save figure
dt = datestr(now, 'yyyymmdd');
print(gcf, img_fmt, ['Fig9_HGP_fXXX_S2_ch', num2str(ch), '_spatialmap', dt])

% Save a color bar image.
figure(); clf
colorbar
colormap(flipud(hot))
caxis(clim)
print(gcf, img_fmt, ['Fig9e_HGP_fXXX_S2_ch', num2str(ch), '_colorbar', dt])

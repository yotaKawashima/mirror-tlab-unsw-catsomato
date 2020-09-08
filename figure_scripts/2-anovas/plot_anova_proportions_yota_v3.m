% Reads data and plots the results of the anova2. 
% Show proportion of channels demonstrating 1) only significant F1 main
% effects. 2) only significant F2 main effects. 3)all in one figure with
% subplot format. 

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%% Setup
% file path
dir_sec = fullfile(data_path, 'collated_data', 'anova_props');
%dir_sec = fullfile(data_path, 'collated_data', ...
%    'anoved_rsampsl_biprref_cmtspwr')

% output image type
%imgtype = '-depsc';
imgtype = '-dsvg';

% Current date and time for output image
dt = datestr(now, 'yyyymmdd');

fsize = 16; % font size

%% Load data to create the matrix of channel proportions

% find filenames
fname = dir(fullfile(dir_sec, ...
    '*_anoved_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_pthresh.mat'));
    
% preallocate
allcat_p_class = [];

% loop across files
for c = 1:numel(fname)
    % load files
    load(fullfile(dir_sec, fname(c).name))
    disp(size(p_class));
    % just truss the new channels onto the end of the old ones
    % dim = (channels x sessions) x frequencies
    % 3 bit (2 main effects and 1 interaction) information is expressed as
    % decimal. e.g. (F2, F1, interaction) = (1,1,1), then 7.
    allcat_p_class = [allcat_p_class; p_class];
end

%% Get proportion data.
% Extract proportions at all frequencies 
prop_all = zeros(7, size(allcat_p_class, 2));
for k = 1:7
    prop_all(k, :) = sum(allcat_p_class == k, 1);
end
prop_all = prop_all / size(allcat_p_class, 1) * 100;

% rearrage order
% 1=001(only interaction), 2=010(only F1), 3=011(F1 and interaction),
% 4=100(only F2), 5=101(F2 and interaction), 6=110(F1 and F2), 7=111(all)
% Then, rearrage into [F1, F2, inter, F1F2, F1inter, F2inter, all].
ytl_order = [2, 4, 1, 6, 3, 5, 7];
prop_all = prop_all(ytl_order, :);
% Add 'Any' class to the proportion matrix. 
% Here, 'Any' class corresponds to all possible combinations of significant
% classes. 
prop_all = [prop_all; sum(prop_all, 1)]; 

% Extract propotions at only frequencies of interest (foi). 
% First get ids of foi.
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);
[~, ind_foi] = parfind_closest(metavars.freq{1}, foi);
% Then, extract propotion at fois.
prop_at_foi = prop_all(:, ind_foi);

%% Plot only F1, only F2 and all with subplot format.
title_list = {'only F1 main effect', 'only F2 main effect', ...
              'only interaction', 'F1 and F2 main effects', ...
              'F1 main effect and interaction', ...
              'F2 main effect and interaction', ...
              'F1 and F2 main effect and interaction', 'any'};

% plot propotion lines
%    1. F1 main interaction
%    2. F2 main interaction
%    3. Interaction
%    4. F1 and F2
%    5. F1 and interaction
%    6. F2 and interaction
%    7. all i.e. F1, F2 and interaction 
%    8. any
% In this script, you only plot 1, 2 and 7. 

figure (1); clf
% Initialisation 
subplotid = 0;

for figid = 1:8
    if figid==1 || figid==2 || figid==7
        subplotid = subplotid + 1;
    else
        continue
    end
    
    subplot(3, 1, subplotid);
    hold on;
    % plot lines (proportion as a function of frequencies)
    % Proportion at all friquences
    plot_val = prop_all(figid,:); 
    q = plot(metavars.freq{1}, plot_val);
    set(q, 'color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    
    % plot propotion as a point only at frequencies of interest.
    % Proportion at frequencies of interest. 
    %{
    plot_val_point = prop_at_foi(figid,:);
    p = scatter(metavars.freq{1}(ind_foi), plot_val_point, 25, 'filled');
    set(p, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    %}
    % plot vertical lines at foi
    
    %   grab axis variables    
    ylims = get(gca, 'YLim');
    ylims = [0, ylims(2)*0.1];
    v1 = numel(foi_f1_and_harm);
    v2 = numel(foi_inter);
    v = zeros(1, numel(foi));
    %   plot vertical lines
    % 23 + harmonics
    v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], ylims, ...
        'Color', [255, 47, 64]/256); 
    % 200
    v(v1+1) = plot([200, 200], ylims, ...
         'Color', [0, 205, 0]/256); 
    % intermodulation
    v(v1+2:v1+v2+1) = plot([foi_inter', foi_inter'], ylims, ...
        'Color', [0, 82, 255]/256); 
    %   format vertical lines
    set(v, 'LineWidth', 2, 'LineStyle', '--');
    hold off

    % format plot area
    ylim auto;
    ylabel('% of significant channels') 
    xlabel('frequency f [Hz]','interpret','none')
    title(title_list(figid));    
    set(gca, 'FontSize', fsize)
    hold off;
end

print(gcf,imgtype,['Fig3_subplot' ,dt])
%% Plot Proportion matrix 
% You can show the matrix by removeing the following curly brackets. 
%{
% Plot proportion of channels per class as a matrix.
figure(2); clf
imagesc(prop_at_foi(1:7,:)) % add sum row

% set labels
xticks = 2:2:numel(foi);
xtickl = foi(xticks);
set(gca, 'xtick', xticks, 'xticklabel', xtickl)
xtickangle(90);

ytl_raw = {'Int', '23', '23 + Int', '200', '200 + Int', '23 + 200', '23 + 200 + Int'};
ytl = ytl_raw(ytl_order);
set(gca, 'yticklabel', ytl)

h = colorbar; 
ylabel(h, '% of significant channels')

xlabel('frequency f [Hz]','interpret','none')
set(gca, 'FontSize', fsize)
print(gcf,imgtype,['Fig3a_' dt])
%}

%% Plot them in one figure.
%{
% plot proportion lines at all frequencies
%   only F1 main effect (corresponding to 2=010)
q(1) = plot(metavars.freq{1}, prop_all(1, :), 'r'); 
%   only F2 main effect (corresponding to 4=100)
q(2) = plot(metavars.freq{1}, prop_all(2, :), 'g');
%   only interaction (corresponding to 1=001)
q(3) = plot(metavars.freq{1}, prop_all(3, :), 'b');
%   sum from 1 to 7, corresponding to the bottom in fig A.
q(4) = plot(metavars.freq{1}, sum(prop_all, 1), 'k');
set(q, 'LineWidth', 1, 'LineStyle', '-')

% plot propotion at fois as points
%   only F1 main effect (corresponding to 2=010)
p(1) = scatter(foi, prop_at_foi(1, :), 25, 'r'); 
%   only F2 main effect (corresponding to 4=100)
p(2) = scatter(foi, prop_at_foi(2, :), 25, 'g'); 
%   only interaction (corresponding to 1=001)
p(3) = scatter(foi, prop_at_foi(3, :), 25, 'b'); 
%   sum from 1 to 7, corresponding to the bottom in fig A. 
p(4) = scatter(foi, sum(prop_at_foi, 1), 25, 'k');
%   format proportion lines
set(p, 'LineWidth', 1)

% plot vertical lines at foi
%   grab axis variables
ylims = get(gca, 'YLim');
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
%   plot vertical lines
% 23 + harmonics
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], ylims, 'r'); 
% 200
v(v1+1) = plot([200, 200], ylims, 'g'); 
% intermodulation
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], ylims, 'b'); 
%   format vertical lines
set(v, 'LineWidth', 1, 'LineStyle', '--')

hold off

% format plot area
ylim([0 20]);
ylabel('% of significant channels') 
xlabel('frequency: f [Hz]','interpret','none')
print(gcf,imgtype,['Fig3b_' dt])
%}
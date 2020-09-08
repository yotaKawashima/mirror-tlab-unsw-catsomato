% Reads data and plots the results of the anova2. 
% Show proportion of channels demonstrating significant 
% main effects and its interaction.
% 1. Show it as a matrix whose row corresponds to a 
% combination of different significant classes. There,
% columns correspond frequencies of interest. 
% (e.g. only F1 main effect, both F1 and F2 main effects,
% etc.)
% 2. Show some rows with continuous frequencies from
% 0Hz to 250Hz.


%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%% Setup
% file path
dir_sec = fullfile(data_path, 'collated_data', 'anova_props');%'anoved_rsampsl_biprref_cmtspwr'

% output image type
%imgtype = '-depsc';
imgtype = '-dsvg';

% Current date and time for output image
dt = datestr(now, 'yyyymmdd');

%% Load data to create the matrix of channel proportions

% find filenames
fname = dir(fullfile(dir_sec, '*_anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_pthresh.mat'));
    
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

%% Plot figure part A
% find the right frequencies and extract frequencies
%{
foi_f1_and_harm = 23:23:250;
[~, ind_f1_and_harm] = parfind_closest(metavars.freq{1}, foi_f1_and_harm);
class_at_foi_f1_and_harm = allcat_p_class(:,ind_f1_and_harm);

foi_f2 = 200;
[~, ind_f2] = parfind_closest(metavars.freq{1}, foi_f2);
class_at_foi_f2 = allcat_p_class(:,ind_f2);

foi_int = sort([177:-23:0 200+23:23:250]);
[~, ind_int] = parfind_closest(metavars.freq{1}, foi_int);
class_at_foi_int = allcat_p_class(:,ind_int);
%}

foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);
[~, ind_foi] = parfind_closest(metavars.freq{1}, foi);

class_at_foi = allcat_p_class(:,ind_foi);

% extract proportions at frequencies of interest.
% We have 3 bit information from i.e. 0~7. But, we do not need 0, which
% expresses that no significance for all 3 (2 main effects and 1
% interaction).
prop_at_foi = zeros(7, size(class_at_foi, 2));
for k = 1:7
    % Sum across (channels*sessions) per frequency.
    prop_at_foi(k, :) = sum(class_at_foi == k, 1);
end
% compute proportion.
prop_at_foi = prop_at_foi / size(class_at_foi, 1) * 100;

% rearrage order
% 2=010(only F1), 4=100(only F2), 1=001(only interaction), 6=110(F1 and
% F2), 3=011(F1 and interaction), 5=101(F2 and interaction), 7=111(all)
% Then, rearrage into [F1, F2, inter, F1F2, F1inter, F2inter, all].
ytl_order = [2, 4, 1, 6, 3, 5, 7];
prop_at_foi = prop_at_foi(ytl_order, :);

colormap jet;
figure(1); clf
imagesc([prop_at_foi; sum(prop_at_foi, 1)]) % add sum row

% set labels
xticks = 2:2:numel(foi);
xtickl = foi(xticks);
set(gca, 'xtick', xticks, 'xticklabel', xtickl)
ytl_raw = {'Int', '23', '23 + Int', '200', '200 + Int', '23 + 200', '23 + 200 + Int'};
ytl = ytl_raw(ytl_order);
ytl{8} = 'any';
set(gca, 'yticklabel', ytl)

h = colorbar; 
title(h, '%')

xlabel('frequency: f [Hz]','interpret','none')
ylabel('Main effects(23,200) and interaction (Int)')

print(gcf,imgtype,['Fig3a_' dt])

%% Plot figure B 

% extract proportions at all frequencies (including frequencies of
% interest)
prop_all = zeros(7, size(allcat_p_class, 2));
for k = 1:7
    prop_all(k, :) = sum(allcat_p_class == k, 1);
end
prop_all = prop_all / size(allcat_p_class, 1) * 100;

% rearrage order
% 2=010(only F1), 4=100(only F2), 1=001(only interaction), 6=110(F1 and
% F2), 3=011(F1 and interaction), 5=101(F2 and interaction), 7=111(all)
% Then, rearrage into [F1, F2, inter, F1F2, F1inter, F2inter, all].
prop_all = prop_all(ytl_order, :);


figure(2); clf
hold on
% plot propotion lines
%    plot 4 lines
%    1. red: 23 and harmonic
%    2. green: 200
%    3. blue: all
%    4. black: any

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
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], ylims, 'r'); % 23 + harmonics
v(v1+1) = plot([200, 200], ylims, 'g'); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], ylims, 'b'); % intermodulation
%   format vertical lines
set(v, 'LineWidth', 1, 'LineStyle', '--')

hold off

% format plot area
ylim([0 20]);
ylabel('% of significant channels') 
xlabel('frequency: f [Hz]','interpret','none')
print(gcf,imgtype,['Fig3b_' dt])

%% Plot the same thing as fig B but for each (F1, F2, inter, and sum)

for fig_id = 1:4
    figure(fig_id + 2); clf
    hold on
    
    switch fig_id
        case 1  % only F1 main effect (corresponding to 2=010)
            plot_val = prop_all(1,:); % at all freq
            plot_val_point = prop_at_foi(1,:); % only at foi
        case 2  % only F2 main effect (corresponding to 4=100)
            plot_val = prop_all(2,:); 
            plot_val_point = prop_at_foi(2,:);
        case 3  % only interaction (corresponding to 1=001)
            plot_val = prop_all(3,:);
            plot_val_point = prop_at_foi(3,:);
        case 4  % sum from 1 to 7, corresponding to the bottom in fig A. 
            plot_val = sum(prop_all,1);
            plot_val_point = sum(prop_at_foi,1);
    end
    % plot lines (proportion as a function of frequencies)
    q = plot(metavars.freq{1}, plot_val);
    set(q, 'LineWidth', 1, 'LineStyle', '-')

    % plot propotion as a point only at frequencies of interest.
    p = scatter(foi, plot_val_point, 25);
    set(p, 'LineWidth', 1);
    
    % plot vertical lines at foi
    
    %   grab axis variables    
    ylims = get(gca, 'YLim');
    v1 = numel(foi_f1_and_harm);
    v2 = numel(foi_inter);
    v = zeros(1, numel(foi));
    %   plot vertical lines
    v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], ylims, 'r'); % 23 + harmonics
    v(v1+1) = plot([200, 200], ylims, 'g'); % 200
    v(v1+2:v1+v2+1) = plot([foi_inter', foi_inter'], ylims, 'b'); % intermodulation
    %   format vertical lines
    set(v, 'LineWidth', 1, 'LineStyle', '--');
    hold off

    % format plot area
    ylim([0 20]);
    ylabel('% of significant channels') 
    xlabel('frequency: f [Hz]','interpret','none')
    print(gcf,imgtype,['Fig3b_', string(fig_id) ,dt])
    title_list = {'only F1 main effect', 'only F2 main effect', 'only interaction', 'sum of any'};
    title(title_list(fig_id));
end

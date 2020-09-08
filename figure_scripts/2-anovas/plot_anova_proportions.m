% Reads data and plots the results of the anova2. 

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
    
    % just truss the new channels onto the end of the old ones
    allcat_p_class = [allcat_p_class; p_class];
end

% find the right frequencies
foi_23 = 23:23:250;
foi_int = [200:-23:0 200+23:23:250];
foi = sort([foi_23 foi_int]);
[~, ind] = parfind_closest(metavars.freq{1}, foi);

% extract frequencies
class_at_foi = allcat_p_class(:, ind);

% extract proportions
prop_at_foi = zeros(7, size(class_at_foi, 2));
for k = 1:7
    prop_at_foi(k, :) = sum(class_at_foi == k, 1);
end
prop_at_foi = prop_at_foi / size(class_at_foi, 1) * 100;

% rearrage order
ytl_order = [2, 4, 1, 6, 3, 5, 7];
prop_at_foi = prop_at_foi(ytl_order, :);


%% Plot figure part A
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

xlabel('response frequency: f [Hz]','interpret','none')
ylabel('2-way ANOVA significance')

print(gcf,imgtype,['Fig3a_' dt])

%% Plot figure part B
figure(2); clf
hold on

% plot propotion lines
%    plot 4 lines
%    1. red: 23
%    2. green: 200
%    3. blue: all
%    4. black: any
p(1) = plot(foi, prop_at_foi(1, :), 'r'); % 23
p(2) = plot(foi, prop_at_foi(2, :), 'g'); % 200
p(3) = plot(foi, prop_at_foi(7, :), 'b'); % int
p(4) = plot(foi, sum(prop_at_foi, 1), 'k');

%   format proportion lines
set(p, 'LineWidth', 2)

% plot vertical lines at foi
%   grab axis variables
ylim = get(gca, 'YLim');
v1 = numel(foi_23);
v2 = numel(foi_int)-1;
v = zeros(1, numel(foi));

%   plot vertical lines
v(1:v1) = plot([foi_23', foi_23'], ylim, 'r'); % 23 + harmonics
v(v1+1) = plot([200, 200], ylim, 'g'); % 200
v(v1+2:v1+v2+1) = plot([foi_int(2:end)', foi_int(2:end)'], ylim, 'b'); % intermodulation

%   format vertical lines
set(v, 'LineWidth', 1, 'LineStyle', '--')

hold off

% format plot area
axis tight
ylabel('% of channel with 2-way ANOVA (q = 0.05)') 
xlabel('response frequency: f [Hz]','interpret','none')

print(gcf,imgtype,['Fig3b_' dt])

%% Plot alternative figure B

% extract proportions at all frequencies
prop_all = zeros(7, size(allcat_p_class, 2));
for k = 1:7
    prop_all(k, :) = sum(allcat_p_class == k, 1);
end
prop_all = prop_all / size(allcat_p_class, 1) * 100;

figure(3); clf
hold on

% plot proportion at all frequencies
%   same format as proportion at foi lines
q(1) = plot(metavars.freq{1}, prop_all(1, :), 'r'); % 23
q(2) = plot(metavars.freq{1}, prop_all(2, :), 'g'); % 200
q(3) = plot(metavars.freq{1}, prop_all(7, :), 'b'); % int
q(4) = plot(metavars.freq{1}, sum(prop_all, 1), 'k');
set(q, 'LineWidth', 0.5, 'LineStyle', '-')

% plot vertical lines at foi
%   grab axis variables
ylim = get(gca, 'YLim');
v1 = numel(foi_23);
v2 = numel(foi_int)-1;
v = zeros(1, numel(foi));

%   plot vertical lines
v(1:v1) = plot([foi_23', foi_23'], ylim, 'r'); % 23 + harmonics
v(v1+1) = plot([200, 200], ylim, 'g'); % 200
v(v1+2:v1+v2+1) = plot([foi_int(2:end)', foi_int(2:end)'], ylim, 'b'); % intermodulation

%   format vertical lines
set(v, 'LineWidth', 1, 'LineStyle', '--')

% plot propotion lines
p(1) = plot(foi, prop_at_foi(1, :), 'r'); % 23
p(2) = plot(foi, prop_at_foi(2, :), 'g'); % 200
p(3) = plot(foi, prop_at_foi(7, :), 'b'); % int
p(4) = plot(foi, sum(prop_at_foi, 1), 'k');

%   format proportion lines
set(p, 'LineWidth', 3 ,'LineStyle', 'none', 'Marker', 'o')


hold off

% format plot area
axis tight
ylabel('% of channel with 2-way ANOVA (q = 0.05)') 
xlabel('response frequency: f [Hz]','interpret','none')

print(gcf,imgtype,['Fig3b_alt_' dt])
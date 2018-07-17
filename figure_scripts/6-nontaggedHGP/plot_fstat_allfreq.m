% Plots the distribution of f-statistic from ANOVA over 0-250 Hz

%% Set overall variables
run(fullfile(mfilename('fullpath'), '../../path_setup.m'))

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';
%data_dir = fullfile(data_path, '/collated_data/anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_fstonly');


%% by area

ylimits = zeros(2, 2);
ax_line = zeros(3, 3);
ax_mark = zeros(3, 2);

ds = 'k-o';

for a = 1:2
    area = ['_S' num2str(a) '_'];
    
    % find names
    cat_names = dirsinside(data_dir);
    for k = 1:numel(cat_names)
        loadname(k) = dir(fullfile(data_dir, cat_names{k}, datatype, ['*' area '*.mat']));
    end
    
    allf = [];
    
    for k = 1:numel(loadname)
        load(fullfile(data_dir, cat_names{k}, datatype, loadname(k).name))
        
        allf = [allf; fstats];
        
    end
    
    % generate mask
    hbw = 0.5;
    f1 = 23;
    f2 = 200;
    fband = [0, 250];
    foi1 = f1:f1:fband(2);
    foi2 = f2:f2:fband(2);
    foi3 = f2:f1:fband(2);
    foi4 = f2:-f1:fband(1);
    foi = sort(unique([foi1 foi2 foi3 foi4]));
    foi = foi(foi>fband(1));
    fs = metavars.freq{1};
    [mask, maskedf] = makefreqmask(fs, foi, fband, hbw);
    
    if a == 1
        mark_ind = floor(linspace(1, length(fs), 23));
    end
    
    % plot
    
    % mean across channels
    fmean = nanmean(allf, 1);
    
    % remove fois
    fmean(1, ~logical(mask), :) = NaN;
    
    % plot
    figure(a); clf;
    hold on
    ax_line(:, a) = plot(fs, fmean(:, :, 1), fs, fmean(:, :, 2), fs, fmean(:, :, 3));
    ax_mark(:, a) = plot(fs(mark_ind), fmean(1, mark_ind, 1), fs(mark_ind), fmean(1, mark_ind, 2), fs(mark_ind), fmean(1, mark_ind, 3));
    hold off
    
    % formatting
    ylimits(:, a) = get(gca, 'YLim');
    xlabel('frequency (Hz)')
    ylabel('f-statistic')
    
end 

offset = 0.2;
ylimfinal(1) = min(ylimits(:, 1))-offset;
ylimfinal(2) = max(ylimits(:, 2))+offset;


%% format area plots
% figure(3); clf;
% ax_line(:, 3) = plot([1, 1], [1, 1], ds, [1, 1], [1, 1], ds, [1, 1], [1, 1], ds);
% legend({'200 Hz effect', '23 Hz effect', 'interaction'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')

% c(3, :) = [180 108 30]/255; % orange
% c(2, :) = [153 204 153]/255; % green
% c(1, :) = [33 131 157]/255; % blue
c = [1,0,0; 0,1,0; 0,0,1];

m(1) = 'o';
m(2) = 'd';
m(3) = 's';
ms = 7;

% set line and marker separately for the areas
for a = 1:2
    figure(a)
    set(gca, 'YLim', ylimfinal)
    
    for k = 1:3
        set(ax_line(k, a), 'Color', c(k, :), 'Marker', 'none', 'LineWidth', 2)
        set(ax_mark(k, a), 'LineStyle', 'none')
    end
    
    print(gcf, '-depsc', ['fig6_S' num2str(a) '_all'])
end

% % set all on line on figure 3
% a = 3;
% 
% figure(a)
% 
% for k = 1:3
%     set(ax_line(k, a), 'Color', c(k, :), 'Marker', m(k), 'MarkerSize', ms)%, 'MarkerEdgeColor', m(k))
% end 
% 
% print(gcf, '-depsc', ['fig6_legend'])
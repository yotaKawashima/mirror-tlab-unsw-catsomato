% Plots the distribution of f-statistic from ANOVA over 0-250 Hz
% Correspond to figure 6 A.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';
%data_dir = fullfile(data_path, '/collated_data/anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_fstonly');

%img_fmt = '-dpng';
img_fmt = '-dsvg';
%% by area

ylimits = zeros(2, 2);
ax_line = zeros(3, 3);
ax_mark = zeros(3, 2);

ds = 'k-o';

for a = 1:2
    area = ['_S' num2str(a) '_'];
    
    % find names
    cat_names = dirsinside(data_dir);
    k=1;
    loadname = dir(fullfile(data_dir, cat_names{k}, datatype, ['*' area '*.mat']));
    for k = 2:numel(cat_names)
        loadname(k) = dir(fullfile(data_dir, cat_names{k}, datatype, ['*' area '*.mat']));
    end
    
    allf = [];
    % load f-stats from all conditions and store it to allf.
    for k = 1:numel(loadname)
        load(fullfile(data_dir, cat_names{k}, datatype, loadname(k).name))
        
        % work out which channels to keep (the bipolar ones)
        if metavars.custom.nsignals == 280 || metavars.custom.nsignals == 176
            lastunipol = prod(metavars.custom.spatialconfig);
            chan = lastunipol+1:metavars.custom.nsignals;
        elseif metavars.custom.nsignals == 180 || metavars.custom.nsignals == 112
            chan = 1:metavars.custom.nsignals;
        else
            error('Only unipolar?');
        end
        
        allf = [allf; fstats(chan, :, :)];
        
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
    
    % plot f-stats (for 200Hz main, 23Hz main, and interaction)
    figure(a); clf;
    hold on
    ax_line(:, a) = plot(fs, fmean(:, :, 1), fs, fmean(:, :, 2), fs, fmean(:, :, 3));
    %ax_mark(:, a) = plot(fs(mark_ind), fmean(1, mark_ind, 1), fs(mark_ind), fmean(1, mark_ind, 2), fs(mark_ind), fmean(1, mark_ind, 3));
    hold off
    
    % formatting
    ylimits(:, a) = get(gca, 'YLim');
    xlabel('frequency [Hz]')
    ylabel('f-statistic [-]')
    
    % compute how many bipolar channels give nan value as f-stats.
    tmp = (1-isnan(allf));
    tmp = tmp(:, 1, 1); 
    fprintf('Number of channels in S%i: %i\n', a, sum(tmp))
    
end 

offset = 0.2;
ylimfinal(1) = min(ylimits(:, 1))-offset;
ylimfinal(2) = max(ylimits(:, 2))+offset;


%% format area plots
% figure(3); clf;
% ax_line(:, 3) = plot([1, 1], [1, 1], ds, [1, 1], [1, 1], ds, [1, 1], [1, 1], ds);
% legend({'200 Hz effect', '23 Hz effect', 'interaction'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')

% (0,1,0) green for F200 main effect
% (1,0,0) red for F23 main effect
% (0,0,1) blue for interaction
c = [0,1,0; 1,0,0; 0,0,1];

%m(1) = 'o';
%m(2) = 'd';
%m(3) = 's';
%ms = 7;

% set line and marker separately for the areas
for a = 1:2
    figure(a)
    set(gca, 'YLim', ylimfinal)
    
    for k = 1:3 % k=1,2,3  (F200, F23, Interaction)
        set(ax_line(k, a), 'Color', c(k, :), 'Marker', 'none', 'LineWidth', 2)
        %set(ax_mark(k, a), 'Color', c(k, :), 'LineStyle', 'none', 'Marker', m(k))
    end
    
    print(gcf, img_fmt, ['Fig6' char(96+a) '_S' num2str(a) '_all'])
end

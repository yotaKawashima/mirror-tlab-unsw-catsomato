% Plots the distribution of f-statistic from ANOVA over 0-250 Hz
% Correspond to figure 6 A.
% Yota modified Rannee's code.

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';
%data_dir = fullfile(data_path, '/collated_data/anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_fstonly');

img_fmt = '-dpng';
%img_fmt = '-dsvg';
%% by area

ylimits = zeros(2, 2);
ax_line = gobjects(3, 2);

ds = 'k-o';
col_l = [[  0, 200,   0];
         [256,   0,   0];
         [  0,   0,   0]]/256; % green, red and blue

col_s = [[170, 256, 170]; 
         [256, 100, 100]; 
         [150, 150, 150]]/256; % green, red and blue
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
    fstd = nanstd(allf, 0, 1);
    
    % remove fois
    fmean(1, ~logical(mask), :) = NaN;
    fstd(1, ~logical(mask), :) = NaN;
    fmean = squeeze(fmean);
    fstd = squeeze(fstd);
    upper = fmean + fstd;
    lower = fmean - fstd; 
    
    % plot f-stats (for 200Hz main, 23Hz main, and interaction)
    figure(a); clf;
    set(gca,'FontSize',16) 
    hold on
    i_fstats = [2,1,3]; % F1main(red), F2main(green), Interact(blue)
    for i_fs = 1:3 % 
        i_fstat = i_fstats(i_fs); % 1 and 2 are flipped. 
        % Interpolate NaN just for showing shade.
        upper_n = upper(:, i_fstat);
        lower_n = lower(:, i_fstat);
        nans = isnan(upper_n);
        upper_n(nans) = interp1(fs(~nans), upper_n(~nans), fs(nans));
        lower_n(nans) = interp1(fs(~nans), lower_n(~nans), fs(nans));
        ax_shaded = fill([fs' fliplr(fs')], ...
                                     [upper_n' fliplr(lower_n')], 'r');
        set(ax_shaded, 'FaceColor', col_s(i_fstat, :));
        set(ax_shaded, 'EdgeColor', col_s(i_fstat, :));
        set(ax_shaded, 'facealpha', .5);
        set(ax_shaded, 'edgealpha', .5);
        ax_line(i_fs, a) = plot(fs, fmean(:, i_fstat), ...
                                'color', col_l(i_fstat,:), ...
                                'Marker', 'none', 'LineWidth', 2);
        
    end    
    %ax_line(:, a) = plot(fs, fmean(:, :, 1), fs, fmean(:, :, 2), fs, fmean(:, :, 3));
    hold off
    
    % formatting
    ylimits(:, a) = get(gca, 'YLim');
    xlabel('frequency [Hz]')
    ylabel('f-statistic [-]')
    
    % compute how many bipolar channels give nan value as f-stats.
    tmp = (1-isnan(allf));
    tmp = tmp(:, 1, 1); 
    fprintf('Number of channels in S%i: %i\n', a, sum(tmp))
    legend(ax_line(:, a), 'F1 main effect', 'F2 main effect', 'Interaction');
end 

offset = 0.2;
%ylimfinal(1) = min(ylimits(:, 1))-offset;
ylimfinal(1) = 0;
ylimfinal(2) = max(ylimits(:, 2))+offset;
for a = 1:2
    figure(a)
    set(gca, 'YLim', ylimfinal)    
    print(gcf, img_fmt, ['Fig6' char(96+a) '_S' num2str(a) '_all'])
end


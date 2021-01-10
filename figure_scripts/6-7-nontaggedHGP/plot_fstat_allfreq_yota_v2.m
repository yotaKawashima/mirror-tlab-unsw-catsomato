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
%img_fmt = '-depsc';
%% by area

ylimits = zeros(2, 2);
ax_line = gobjects(3, 2);

ds = 'k-o';     
     
% Store subplot axis for lim
hAxs = cell(2,1);
titles = {'F1 main effect', 'F2 main effect', 'Interaction'};

% Store mean f-stats for each area
fmeans = cell(2,1);

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
    % Remove 50, 100 and 150Hz as well. (Remove line noise freq.)
    foi = sort(unique([foi1 50 100 150 foi2 foi3 foi4]));
    foi = foi(foi>fband(1));    
    fs = metavars.freq{1};
    [mask, maskedf] = makefreqmask(fs, foi, fband, hbw);
   
    % plot
    
    % mean across channels
    fmean = nanmean(allf, 1);
    %fstd = nanstd(allf, 0, 1);
    quantile1 = prctile(allf, 25, 1);
    quantile2 = prctile(allf, 50, 1);
    quantile3 = prctile(allf, 75, 1);
    prct1 = prctile(allf, 90, 1); 
    
    % remove fois
    fmean(1, ~logical(mask), :) = NaN;
    %fstd(1, ~logical(mask), :) = NaN;
    quantile1(1, ~logical(mask), :) = NaN;
    quantile2(1, ~logical(mask), :) = NaN;
    quantile3(1, ~logical(mask), :) = NaN;
    prct1(1, ~logical(mask), :) = NaN;
  
    fmean = squeeze(fmean);
    fmeans{a, 1} = fmean; % Store data for histogram
    %fstd = squeeze(fstd);
    quantile1 = squeeze(quantile1);
    quantile2 = squeeze(quantile2);
    quantile3 = squeeze(quantile3);
    prct1 = squeeze(prct1);
    
    
    % plot f-stats (for 200Hz main, 23Hz main, and interaction)
    figure(a); clf;
    set(gca,'FontSize',16) 
    %Preallocate subplots
    hAx = gobjects(3,1);
    
    i_fstats = [2,1,3]; % F1main(red), F2main(green), Interact(blue)

    for i_fs = 1:3 % 
        hAx(i_fs) = subplot(1, 3, i_fs);
        hold on
        i_fstat = i_fstats(i_fs); % 1 and 2 are flipped. 

        % Plot quantile
        plot(fs, quantile1(:, i_fstat), 'color', [200, 200, 200]/256,  ...
             'Marker', 'none', 'LineWidth', 1);
        % Plot mean
        plot(fs, fmean(:, i_fstat), 'color', 'k', ...
             'Marker', 'none', 'LineWidth', 2);

        % Median
        plot(fs, quantile2(:, i_fstat), 'color', [200, 200, 200]/256, ...
            'Marker', 'none', 'LineWidth', 1);
         
        % Quantile
        plot(fs, quantile3(:, i_fstat), 'color', [200, 200, 200]/256, ...
            'Marker', 'none', 'LineWidth', 1);

        plot(fs, prct1(:, i_fstat), 'color', [200, 200, 200]/256, ...
            'Marker', 'none', 'LineWidth', 1);
        
        hold off
        
        legend('25%', 'Mean', '50%', '75%', '90%');
        title(titles{i_fs});
        xlabel('frequency [Hz]')
        ylabel('f-statistic [-]')
    
    end
    
    % formatting
    max_y = max(nanmax(prct1));
    min_y = min(nanmin(quantile1));
    ylimits(:, a) = [min_y, max_y];
    
    % compute how many bipolar channels give nan value as f-stats.
    tmp = (1-isnan(allf));
    tmp = tmp(:, 1, 1); 
    fprintf('Number of channels in S%i: %i\n', a, sum(tmp))
    
    % Store subplot axis
    hAxs{a,1} = hAx;
    
end 

offset = 0.2;
%ylimfinal(1) = min(ylimits(:, 1))-offset;
ylimfinal(1) = 0;
ylimfinal(2) = max(ylimits(:, 2))+offset;
for a = 1:2
    figure(a)
    hAx = hAxs{a,1};
    set(hAx, 'YLim', ylimfinal);  
    set(hAx, 'Xlim', [0 250]);
    %print(gcf, img_fmt, ['Fig6' char(96+a) '_S' num2str(a) '_all'])
end


% Plot histogram f-stats at 100Hz
[~, f_ind] = find_closest(metavars.freq{1,1}, 100);


for a = 1:2
    figure(2+a);clf
    hold on
    for i_fs = 1:3 % 
        i_fstat = i_fstats(i_fs); % 1 and 2 are flipped. F1, F2, interaction
        fmean = fmeans{a, 1};
        h = histogram(log(fmean(:, i_fstat)+1));%/sum(fmean(:, i_fstat)))
        h.BinWidth = 0.02;
    end
    xlabel('f-stats at 100Hz');
    ylabel('# of data points');
    title(['S', num2str(a), ' histogram']);
    legend('F1 main effect', 'F2 main effect', 'Interaction');
    hold off
    print(gcf, img_fmt, ['Fig6_S' num2str(a) '_histogramAt100HZ'])

end


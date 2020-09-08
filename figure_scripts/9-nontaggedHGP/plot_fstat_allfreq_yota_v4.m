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

%img_fmt = '-dpng';
img_fmt = '-depsc';
%% by area
ylimits = zeros(2, 2);
% Store subplot axis for lim
hAxs = cell(2,1);
titles = {'F1 main effect', 'F2 main effect', 'Interaction'};

% Store mean f-stats for each area
fmeans = cell(2,1);

% Plot S1 and S2 in the same figure
%figure(1);clf
colors = cell(2,1);
colors{1,1} = [  0,   0,   0; ...
               125, 125,  125; ...
               220, 220, 220]/ 256;
colors{2,1} = [255,   0,   0; ...
               255,  90,  90; ...
               255, 180, 180]/ 256;
 
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
    sessions = [];
    bpchannel_matrix = [];
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
        
        allf = [allf; fstats(chan, :, :)]; % channels x freqs x 3
        sessions = [sessions ; ones(length(chan),1)*k];
        bpchannel_matrix = [bpchannel_matrix, 1:length(chan)];
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
    prct1 = prctile(allf, 95, 1);
    prct2 = prctile(allf, 90, 1);
    prct3 = prctile(allf, 50, 1);
    
    % remove fois
    prct1(1, ~logical(mask), :) = NaN;
    prct1 = squeeze(prct1);
    prct2(1, ~logical(mask), :) = NaN;
    prct2 = squeeze(prct2);
    prct3(1, ~logical(mask), :) = NaN;
    prct3 = squeeze(prct3);
       
    % plot f-stats (for 200Hz main, 23Hz main, and interaction) after 
    % smooting. (x-axis, y-axis) = (frequency, f-stats)
    % Smooth plots by removing fluctuation within 1Hz width. (This does not
    % mean we remove 1Hz wave components of lines.)
    i_fstats = [2,1,3]; % F1main, F2main, Interact
    % Get # of samples within 1Hz. (i.e. samples/1Hz). 
    [~, samples_per_1Hz] = find_closest(fs, 1);
    
    for i_fs = 1:3 % 
        figure(i_fs);
        %hAx(i_fs) = subplot(1, 3, i_fs);
        hAx(i_fs) = gca;
        hold on
        i_fstat = i_fstats(i_fs); % 1 and 2 are flipped. 
        
        % Smoothing for visualisation.
        % First, interpolate NaN for smoothing.
        prct1_n = prct1(:, i_fstat);
        prct2_n = prct2(:, i_fstat);
        prct3_n = prct3(:, i_fstat);
        nans = isnan(prct1_n);
        prct1_n(nans) = interp1(fs(~nans), prct1_n(~nans), fs(nans));
        prct2_n(nans) = interp1(fs(~nans), prct2_n(~nans), fs(nans));
        prct3_n(nans) = interp1(fs(~nans), prct3_n(~nans), fs(nans));
   
        % Second, filter out components higher than (1sample/1Hz)
        prct1_lpass = lowpass(prct1_n, 1, samples_per_1Hz);
        prct2_lpass = lowpass(prct2_n, 1, samples_per_1Hz);
        prct3_lpass = lowpass(prct3_n, 1, samples_per_1Hz);
        
        % Smoothed line
        colors_now = colors{a,1};
        % Plot percentile
        plot(fs, prct1_lpass, 'color', colors_now(1,:),  ...
             'Marker', 'none', 'LineWidth', 2);
        % Median
        plot(fs, prct2_lpass, 'color', colors_now(2,:), ...
            'Marker', 'none', 'LineWidth', 2);
        % Quantile
        plot(fs, prct3_lpass, 'color', colors_now(3,:), ...
            'Marker', 'none', 'LineWidth', 2);
        
        hold off
        
        if a == 2
            legend('S1:Top 5%', 'S1:Top 10%', 'S1:Top 50%', ...
                   'S2:Top 5%', 'S2:Top 10%', 'S2:Top 50%');
          
            title(titles{i_fs});
            xlabel('frequency [Hz]')
            ylabel('f-statistic [-]')
        end %if a==2
    end
    
    % Get limits for formatting
    max_y = max(nanmax(prct1));
    min_y = min(nanmin(prct3));
    ylimits(:, a) = [min_y, max_y];
    
    % compute how many bipolar channels give nan value as f-stats.
    % Need to check only one freq and one f-stat for each channel.
    invalid_chs = isnan(allf(:, 1, 1));
    num_invalid_chs = sum(invalid_chs);
    if num_invalid_chs == 0
        fprintf('Number of valid channels in S%i: %i out of %i\n', ... 
                a, size(allf,1) - num_invalid_chs ,size(allf,1))
        fprintf('No invalid bipolar chs\n');
    else % if there are invalid bipolar channels
        fprintf('Number of valid channels in S%i: %i out of %i\n', ... 
                a, size(allf,1) - num_invalid_chs ,size(allf,1))

        % Get number of invalid bipolar channels for each sessoin.
        sessionid_invalidch = sessions(invalid_chs);
        bpchannel_invalidch = bpchannel_matrix(invalid_chs);
        unique_sessionids = unique(sessionid_invalidch);  
        for i_session = 1:length(unique_sessionids)
            session_id = unique_sessionids(i_session);
            fprintf('%s -- # of invalid bipolar chs: %i\n', ...
                    loadname(session_id).name(1:13),...
                    sum(sessionid_invalidch==session_id));
            disp('bipolar channel id: ');
            disp(bpchannel_invalidch((sessionid_invalidch == session_id)));
        end % i_session = 1:length(unique_sessionids)
    end % num_invalid_chs == 0
end 

offset = 0.2;
%ylimfinal(1) = min(ylimits(:, 1))-offset;
ylimfinal(1) = 0;
ylimfinal(2) = max(ylimits(:, 2))+offset;

for i_fs = 1:3
    figure(i_fs);
    set(gca, 'YLim', ylimfinal);  %set(hAx, 'YLim', ylimfinal);  
    set(gca, 'Xlim', [0 250])     %set(hAx, 'Xlim', [0 250]);
    set(gca,'FontSize',16);    %set(hAx,'FontSize',16);
    print(gcf, img_fmt, ['Fig9_' titles{i_fs} '_all'])
end % i_fs = 1:3

% Create color bar as legend.
for a=1:2
    figure();clf
    ax = gca;
    ax.Visible = 'off';
    colormap(flipud(colors{a,1}));
    c = colorbar;
    caxis([0, 1]);
    c.Location = 'west';
    c.Ticks = [1/6, 0.5, 1-1/6];
    c.TickLabels = {'top 50%', 'top 10%', 'top 5%'};
    set(gca,'FontSize',16);    %set(hAx,'FontSize',16);
    print(gcf, img_fmt, ['Fig9_S' num2str(a) '_all_colorbar']) 
end % a=1:2
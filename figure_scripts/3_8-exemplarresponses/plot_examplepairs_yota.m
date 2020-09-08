% Plots example response across conditons from Session 2-1. 
% You can plot expemplar responese of three different types.
% Plot both log SNR and VELogP.

%% Data selection
cat_name = 'C20110808_R03';

% area (S1 or S1), channel, frequency
plotpair{1} = [2, 49, 200]; 
plotpair{2} = [1, 51, 23];
plotpair{3} = [1, 33, 200];

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%
%img_fmt = '-dpng';
img_fmt = '-dsvg';

%% Set script specific variables
%logpower_file = 'epoched_rsampsl_biprref_evkresp_cmtspwr';
logsnr_file = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';
velogp_file = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr';

%% Find files
% find list of names
data_dir = fullfile(data_path, 'included_datasets');
loadnames = dir(fullfile(data_dir, cat_name, logsnr_file, '*.mat'));

% find some numbers
nCond = numel(loadnames)/2; % number of conditions per area 
allps = cell(1, 2); % one for each area
allhs = cell(1, 2);
meansnrs = cell(1, 2);

for exampleid=1:length(plotpair)
    area = plotpair{1,exampleid}(1);
    bipolar_channel = plotpair{1,exampleid}(2);
    foi = plotpair{1,exampleid}(3);
    % Loop through conditions in the session.
    for condid = 1:nCond
        
        % Load logsnr data.
        fprintf('Loading area %i/%i, loading data %2i/%i\n', ...
            area, 2, condid, nCond)
        load(fullfile(data_dir, cat_name, logsnr_file, ...
            loadnames((area-1)*nCond + condid).name))
        
        % if it's the first condition, 
        % 1) check whether data included unipolar, and
        % 2) preallocate frequencies of interest (per area).
        if condid==1
            % 1) check if unipolar id included in srn data. If included,
            %   unipolar should be removed. 
            if strcmp(data.label{1}(1:3), 'raw')
                % the first bipolar channel
                bipolarchid = 1 + prod(data.custom.spatialconfig);
                % # of bipolar channels
                nChan = data.custom.nsignals - bipolarchid + 1; 
            else
                nChan = data.custom.nsignals;
                bipolarchid = 1;
            end
            
            % 2) preallocate data
            allps = zeros(nChan, nCond); % bipolar ch x conditions
            allhs = zeros(nChan, nCond);
            means_forspatialmap = zeros(nChan, nCond);
    
            % Get id and actual value of foi in continuous data
            [foival, foiid] = find_closest(data.freq{1}, foi);
            
        else % condid>1
            % check that the value of frequency is right
            if ~isequal(data.freq{1}(foiid), foival)
                warning(['Frequency indicies inconsistent.', ... 
                    ' Recalculating. (a=%i, c=%i)'], area, condid)
    
                [foiid, foival] = find_closest(data.freq{1}, foi);
            end
        end % if c=1
        
        % Extract log SNR data around foi from all bipolar channels.
        % x axis data for single channel plot 
        % frequency constants
        [~, f_i] = find_closest(data.freq{1}, foi);
        iss = mean(diff(data.freq{1}));
        hbw_inds = ceil(3/iss);
        freqs = data.freq{1}(f_i-hbw_inds:f_i+hbw_inds);
        
        % Initialise storage for single channel
        if condid ==1
            means_forsinglech = zeros(length(freqs), nCond);
            stds_forsinglech = zeros(length(freqs), nCond);
            freqs_forsinglech = zeros(length(freqs), nCond);
        end
        
        % Store freqs for each condition just for plotting.
        % freqs should be the same across conditions. 
        freqs_forsinglech(:, condid) = freqs;
        
        logsnrs_aroundfoi = squeeze(...
                data.trial(bipolar_channel, f_i-hbw_inds:f_i+hbw_inds, :));
        
        % mean and std at a bipolar channel across trials.
        means_forsinglech(:, condid) = mean(logsnrs_aroundfoi, 2);
        stds_forsinglech(:, condid) = std(logsnrs_aroundfoi, 0, 2);
        
        % Extract log SNR data at foi from all bipolar channels.
        logsnrs_atfoi = squeeze(data.trial(bipolarchid:end, foiid, :)); 
        % mean across trials within a condition.
        means_forspatialmap(:, condid) = mean(logsnrs_atfoi, 2);
        
        % t-test per channel. Variance: trials within a stimulus condition. 
        for ch = 1:nChan % channel
            [allhs(ch,condid) , allps(ch,condid)] = ...
                ttest(logsnrs_atfoi(ch, :), 0, 'Tail', 'right');
        end %ch=1:nChan
               
    end % condid=1:nCond

    % Threshold SNR
    q = 0.05;

    %means(allhs==0)=NaN;

    % pID is a threshold for p-val. (scalar not vector)
    [pID, ~] = eeglab_fdr(allps(:, :), q, 'parametric');

    % Overwrite mean data: set NaN for not significant channels.
    tmp = means_forspatialmap; % channel x conditions.
    tmp(allps(:, :)>=pID)=NaN; 
    means_forspatialmap = tmp;
    fprintf('S%i f=%3i: %f\n', area, foiid, pID)

    labels = data.label(bipolarchid:end); % label of bipolar channels
    spatialconfig = data.custom.spatialconfig; % channel config
    condfig = data.custom.subplotconfig; % condition config
    
    params_filename = {area, bipolar_channel, foi, condfig, img_fmt};    
   
    % Plot single channel data
    plotsinglechannels(freqs_forsinglech, ...
        means_forsinglech, stds_forsinglech, params_filename);
   
    % Check whether means_forspatialmap is the corresponds to single
    % channel data.
    figure();
    scatter(linespace(1, size(means_forespatialmap, 2)), means_forespatialmap(bipolar_channel,:));
    pause;
    clf;
    
    % Plot spatial map 
    plotspatialmap(means_forspatialmap, spatialconfig, ...
        labels, params_filename);
  
end % exampleid

%% Plot single channels
function plotsinglechannels(x_data_array, means, stds, params_filename)
    % Plot single channel data per stimulus condition.
    % Input 
    %   x_datas           : data for x axis
    %   means             : mean data 
    %   stds              : std data 
    %   params_filename   : parameters for file name 
    %           area            : S1 or S2 
    %           bipolar_channel : bipolar channel id
    %           foi             : frequency of interest
    %           condfig         : stimulus condition 
    %           img_fmt         : image format
    %   
    % Output 
    %   figure 
    figure(); clf
    
    area = params_filename{1};
    bipolar_channel = params_filename{2};
    foi = params_filename{3};
    condfig = params_filename{4};
    img_fmt = params_filename{5};
    
    nConds = condfig(1) * condfig(2);
    
    % color for shaiding
    col = [192,192,192]/256;
    
    % For shading
    lowers = means - stds;
    uppers = means + stds;
    
    % Preallocate subplots.
    hAx = gobjects(4*4,1);
    
    for condid = 1:nConds
        % Subplot.
        hAx(condid) = subplot(4, 4, condid);
        
        % Shading 
        lower = lowers(:, condid)';
        upper = uppers(:, condid)';
        mean = means(:, condid)';
        x_data = x_data_array(:, condid)';
        hold on 
        h1 = fill([x_data fliplr(x_data)], [upper fliplr(lower)], 'r');
        %h2 = fill([x_data flx_data(end) x_data(1)], [lower 0 0], 'r');
        set(h1, 'FaceColor', col)
        set(h1, 'EdgeColor', col)
        %set(h2, 'FaceColor', col)
        %set(h2, 'EdgeColor', col)
        
        % Line
        plot(x_data, means(:, condid), 'color', 'k', ...
            'LineWidth', 1);
        hold off
 
    end
    
    % Get min and max across stimulus conditions 
    min_ = min(min(lowers));
    max_ = max(max(uppers));
    lims = [min_-1 max_+1];
    
    % Set x and y axis
    xticks(hAx, [foi-3 foi foi+3]);
    xlim(hAx, [foi-3 foi+3]);
    ylim(hAx, lims);

    % save figure
    print(gcf, img_fmt, ['Fig4-5_f' num2str(foi) '_S' num2str(area) ...
                         '_ch' num2str(bipolar_channel)])

end

%% Plot spatial map
function plotspatialmap(data, spatialconfig, labels, params_filename)
    % Plot single channel data per stimulus condition.
    % Input 
    %   data              : data shown as color
    %   spatial config    : spatial configuration 
    %   labels            : labels 
    %   params_filename   : parameters for file name 
    %           area            : S1 or S2 
    %           bipolar_channel : bipolar channel id
    %           foi             : frequency of interest
    %           condfig         : stimulus condition 
    %           img_fmt         : image format
    %   
    % Output 
    %   figure 
    figure(); clf
    
    area = params_filename{1};
    bipolar_channel = params_filename{2};
    foi = params_filename{3};
    condfig = params_filename{4};
    img_fmt = params_filename{5};
    
    nConds = condfig(1) * condfig(2);
    
    % Relabel for the draw_biprref function
    relabels = draw_biprref_chlabfunction(labels);
    
    % find channel of interest location
    [x1, y1] = patch_helper(bipolar_channel, relabels);

    % Set colormap
    %clim = find_lims(data);
    %cmap = [[240 240 240]/256; flipud(cool(ceil(diff(clim))*10))];
    %cmap = [[240 240 240]/256; cool(ceil(diff(clim))*10)];
    %colormap default;
    colormap(flipud(hot));

    % Preallocate subplots.
    hAx = gobjects(4*4,1);
    
    % Remove inf and cet color lim.
    data_ = data(~isinf(data));
    clim = [nanmin(data_) nanmax(data_)];

    for condid = 1:nConds
        % make subplot
        hAx(condid) = subtightplot(condfig(1), condfig(2), condid);
        
        % goes across rows
        values = data(:, condid);
        %[p, FV] = draw_biprref_vy(values, relabels, spatialconfig);
        draw_biprref_vy(values, relabels, spatialconfig, clim);
        hold on
        % add patch to foi
        p1 = patch(x1, y1, 12);
        set(p1, 'EdgeColor', [0, 0, 200]/256, 'LineWidth', 1.5, 'FaceColor', 'none')
        hold off
    end
    
    % Set the same clim across subplots.
    %set(hAx, 'CLim', clim)
    
    % Make axsis invisible
    set(hAx, 'XTicklabel',[])
    set(hAx, 'YTicklabel',[])
    
    % save figure
    print(gcf, img_fmt, ['Fig4-5_f' num2str(foi) '_S' num2str(area)])
    
    % Save a color bar image.
    figure(); clf
    colorbar
    colormap(flipud(hot))
    caxis(clim)
    print(gcf, img_fmt, ...
        ['Fig4-5_f' num2str(foi) '_S' num2str(area) '_colorbar'])
end
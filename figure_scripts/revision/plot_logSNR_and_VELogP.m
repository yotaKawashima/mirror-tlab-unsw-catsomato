% Plot logSNR and VELogP to check that logSNR is smoothed by its defintion
% (i.e. subtracting responses of neighboring freqs from the one at frequency 
% of our interest) compared to VELogP
% Focus on one channel in the maximum stimulus amplitude condition.

%% Data selection
cat_name = 'C20110808_R03';

% area (S1 or S1), channel, frequency, type(1:logSNR, 2:VELogP)
% Note : pval from ANOVA (F1 main, F2 main, interaction) 
% logSNR
% Fig 3: only F1 main effect
plotpair{1} = [1, 131, 23,  1]; %(0.000000, 0.053872, 0.517583)
% Fig 4: only F2 main effect
%plotpair{2} = [1, 36, 200,  1]; %(0.000425, 0.000000, 0.134479)
% Fig 5: all 
%plotpair{3} = [1, 158, 23,  1]; %(0.000000, 0.000000, 0.000000)

% VELogP
% Fig 6: only F1 main effect
plotpair{2} = [1, 131,  23, 2]; %(0.000000, 0.195915, 0.041284)
% Fig 7: only F2 main effect
%plotpair{5} = [1, 122, 200, 2]; %(0.003258, 0.000000, 0.044181) 
% Fig 8: all 
%plotpair{6} = [1, 43,  23,  2]; %(0.000000, 0.000000, 0.000000) 


%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));
%
%img_fmt = '-depsc2';
%img_fmt = '-dpng';
%img_fmt = '-dtiff';
img_fmt = '-dpdf';
fsize = 10; % font size
ftype = 'Arial'; % font type
x_width = 8; % fig width
y_width = 8; % fig height

%% Set script specific variables
%logpower_file = 'epoched_rsampsl_biprref_evkresp_cmtspwr';
dtypes = {'snrsurr', 'evkdpwr'};

%% Find files
% find list of names
data_dir = fullfile(data_path, 'included_datasets');

% figure 
figure(); clf
% Preallocate subplots.
hAx = gobjects(2,1);

for exampleid=1:length(plotpair)
    area = plotpair{1,exampleid}(1);
    bipolar_channel = plotpair{1,exampleid}(2);
    foi = plotpair{1,exampleid}(3);
    dattype = dtypes{plotpair{1,exampleid}(4)};
    file_type = ['epoched_rsampsl_biprref_evkresp_cmtspwr_', dattype];
    
    % Find files
    loadnames = dir(fullfile(data_dir, cat_name, file_type, '*.mat'));
    % Get the max stimulus vibration condition.
    nCond = numel(loadnames)/2; % number of conditions per area 
        
    % Load logSNR/VELogP data.
    fprintf('Loading area %i, loading data %s\n', area, dattype)
    load(fullfile(data_dir, cat_name, file_type, ...
        loadnames((area-1)*nCond + nCond).name))

    % 1) check whether data included unipolar
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

    % 2) preallocate frequencies of interest (per area).
    %allps = zeros(nChan, nCond); % bipolar ch x conditions
    %allhs = zeros(nChan, nCond);
    means_forspatialmap = zeros(nChan, nCond);

    % Get id and actual value of foi in continuous data
    [foival, foiid] = find_closest(data.freq{1}, foi);

    % Extract log SNR data around foi from all bipolar channels.
    % x axis data for single channel plot 
    % frequency constants
    iss = mean(diff(data.freq{1}));
    hbw_inds = ceil(3/iss);
    freqs = data.freq{1}(foiid-hbw_inds:foiid+hbw_inds);
        
    % Initialise storage for single channel
    if exampleid ==1
        means_forsinglech = zeros(length(freqs), 2);
        stds_forsinglech = zeros(length(freqs), 2);
        freqs_forsinglech = zeros(length(freqs), 2);
    end

    % Store freqs for each condition just for plotting.
    % freqs should be the same across conditions. 
    freqs_forsinglech(:, exampleid) = freqs;

    % data.trial = channels x frequencies x trials 
    data_bipolar = data.trial(bipolarchid:end,:,:);
    % mean and std at bipolar channels and frequencies across trials.
    data_bipolar_mean = mean(data_bipolar, 3);
    data_bipolar_std = std(data_bipolar, 0, 3);

    % For single channel around foi
    means_forsinglech(:, exampleid) = data_bipolar_mean(....
        bipolar_channel, foiid-hbw_inds:foiid+hbw_inds);
    stds_forsinglech(:, exampleid) = data_bipolar_std(...
        bipolar_channel, foiid-hbw_inds:foiid+hbw_inds);

    % For all channels at foi
    means_forspatialmap(:, exampleid) = data_bipolar_mean(:,foiid);
    
    labels = data.label(bipolarchid:end); % label of bipolar channels
    
    % color for shaiding
    col = [192,192,192]/256;
    
    % For shading
    lowers = means_forsinglech - stds_forsinglech;
    uppers = means_forsinglech + stds_forsinglech;

    
    % Subplot.
    hAx(exampleid) = subplot(2, 1, exampleid);

    % Shading 
    lower = lowers(:, exampleid)';
    upper = uppers(:, exampleid)';
    x_data = freqs_forsinglech(:, exampleid)';
    hold on 
    h1 = fill([x_data fliplr(x_data)], [upper fliplr(lower)], 'r');
    %h2 = fill([x_data flx_data(end) x_data(1)], [lower 0 0], 'r');
    set(h1, 'FaceColor', col)
    set(h1, 'EdgeColor', col)
    %set(h2, 'FaceColor', col)
    %set(h2, 'EdgeColor', col)

    % Line
    plot(x_data, means_forsinglech(:, exampleid), 'color', 'k', ...
        'LineWidth', 1);
    % Set y label
    if plotpair{1,exampleid}(4) == 1
        ylabel_text = 'logSNR [dB]';
        % Add line showing 0.5Hz range 
        plot([foi-0.25, foi+0.25], [7, 7], ...
            'color', 'k',...
            'linewidth', 1.5);
        text(24, 7, '0.5Hz width')
        % Set x axis
        ylim([-5, 8]);
    elseif plotpair{1,exampleid}(4) == 2
        ylabel_text = 'VELogP [dB]';
        % Add line showing 0.5Hz range 
        plot([foi-0.25, foi+0.25], [17, 17], ...
            'color', 'k',...
            'linewidth', 1.5);
        text(24, 17, '0.5Hz width')

        % Set x axis
        ylim([-15, 19]);
    end
    ylabel(ylabel_text)
    grid on
    grid minor
    hold off
   
end % exampleid
%%
% Set x axis
xticks(hAx, foi-3:1:foi+3);
xlim(hAx, [foi-3, foi+3]);
xlabel(hAx, 'frequency [Hz]');

% save figure
filename = ['FigRev4_logSNR_VELogP_f' num2str(foi) '_S' ...
            num2str(area) '_ch' num2str(bipolar_channel)];
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');    
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, filename);
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, filename, '-r300');
end

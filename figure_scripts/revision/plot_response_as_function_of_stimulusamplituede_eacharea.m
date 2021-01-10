% Plot logSNR at foi as a function of stimulation amplitude.
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
x_width = 25; % fig width
y_width = 25; % fig height

% frequencis of our interest
% (fundamental, harmonics, intermodulations)
f1 = 23;
f2 = 200;
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
%foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);
foi = [f1, f2];

% Data directory
data_dir = fullfile(data_path, 'included_datasets');

%% First focus on only one channel in a session
%{
cat_name = 'C20110808_R03';
area = 1;
file_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'; %evkdpwr
channel_id = 131; % channel we plot

% Find files
loadnames = dir(fullfile(data_dir, cat_name, file_type, '*.mat'));
nCond = numel(loadnames)/2; % number of conditions per area 

% Loop through conditions in the session.
for condid = 1:nCond

    % Load logSNR/VELogP data.
    fprintf('Loading area %i/%i, loading data %2i/%i\n', ...
        area, 2, condid, nCond)
    load(fullfile(data_dir, cat_name, file_type, ...
        loadnames((area-1)*nCond + condid).name))
    
    % if it's the first condition, 
    % 1) check whether data included unipolar, and
    % 2) preallocate frequencies of interest (per area).
    % 3) preallocate stimulus amplitude.
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
        nTrials = size(data.trial, 3);
        resp_matrix = zeros(nChan, length(foi), ...
                            data.custom.subplotconfig(1),...
                            data.custom.subplotconfig(2),...
                            nTrials);
        
        % 3) preallocate stim info matrix
        stim_f1_amplitudes = zeros(size(data.custom.subplotconfig));
        stim_f2_amplitudes = zeros(size(data.custom.subplotconfig));
        % Get id and actual value of foi in continuous data
        [foivals, foiids] = find_closest(data.freq{1}, foi);

    else % condid>1
        % check that the value of frequency is right
        if ~isequal(data.freq{1}(foiids(1)), foivals(1))
            warning(['Frequency indicies inconsistent.', ... 
                ' Recalculating. (a=%i, c=%i)'], area, condid)

            [foivals, foiids] = find_closest(data.freq{1}, foi);
        end
    end % if c=1
    
    % Set amplitude id
    f1_id = ceil(condid / data.custom.subplotconfig(2));
    f2_id = mod(condid, data.custom.subplotconfig(2));
    if f2_id == 0
        f2_id = data.custom.subplotconfig(2);
    end % f2_ind == 0
    % get amplitude information from file name
    stim_f1_amplitudes(f1_id, f2_id) = str2double(loadnames((area-1)*nCond + condid).name(23:25));
    stim_f2_amplitudes(f1_id, f2_id) = str2double(loadnames((area-1)*nCond + condid).name(32:34));
    
    % data.trial = channels x frequencies x trials 
    data_bipolar = data.trial(bipolarchid:end,:,:);
    
    % Extract responses at foi (mean across trials)
    % channels x freqs x F1 stim amp x F2 stim amp x trials
    resp_matrix(:, :, f1_id, f2_id , :) = data_bipolar(:, foiids, :);

end % condid=1:nCond

% mean/std resp at f1 and f2 across trials at the specified bipolar channel
mean_respf1_matrix = ...
    squeeze(mean(resp_matrix(channel_id, (foi==f1), :, :, :), 5));
stderr_respf1_matrix = ...
    squeeze(std(resp_matrix(channel_id, (foi==f1), :, :, :), 0, 5)) ...
    / sqrt(nTrials);
mean_respf2_matrix = ...
    squeeze(mean(resp_matrix(channel_id, (foi==f2), :, :, :), 5));
stderr_respf2_matrix = ...
    squeeze(std(resp_matrix(channel_id, (foi==f2), :, :, :), 0, 5)) ...
    / sqrt(nTrials);
% Plot responses as a function of stimulus amplitudes
figure();


% Legend label
n_f1_amp = size(stim_f1_amplitudes, 1);
n_f2_amp = size(stim_f1_amplitudes, 2);
legend_label_f1 = {};
legend_label_f2 = {};
for f_i = 1:n_f1_amp
    legend_label_f1{f_i} = [num2str(stim_f1_amplitudes(f_i, 1)), '\mum'];
    legend_label_f2{f_i} = [num2str(stim_f2_amplitudes(1, f_i)), '\mum'];
end

% Set coloar
cmp_tmp = colormap('parula');%gray');
%color_step = floor((size(cmp_tmp, 1) - 40)/n_f1_amp);
%cmp = flipud(cmp_tmp(1:color_step:256-40, :));
color_step = floor(size(cmp_tmp-10, 1)/n_f1_amp);
cmp = flipud(cmp_tmp(1:color_step:256-10, :));

% plot 23Hz responses
hAx1 = subplot(1, 2, 1);
hold on

% Use different color for different stim f2 amps.
h_offset = [-0.3, -0.15, 0.15, 0.3];
for f2_amp_id = 1:n_f2_amp    
    x_values = 1:n_f1_amp;
    %x_values = x_values + h_offset(f2_amp_id);
    errb = errorbar(x_values, mean_respf1_matrix(:, f2_amp_id), ...
             stderr_respf1_matrix(:, f2_amp_id), ...
             'color', cmp(f2_amp_id, :));
    % Remove legend for error bar
    errb_annotation = get(errb, 'Annotation');
    set(get(errb_annotation,'LegendInformation'),...
            'IconDisplayStyle','off');

    scatt = scatter(x_values, mean_respf1_matrix(:, f2_amp_id), ...
                    6^2, cmp(f2_amp_id, :), 'filled');
    % Remove legend for scatter
    %scatt_annotation = get(scatt, 'Annotation');
    %set(get(scatt_annotation,'LegendInformation'),...
    %        'IconDisplayStyle','off');
end
hold off
xlim([0 n_f1_amp+1])
xticks(1:n_f1_amp)
xticklabels(stim_f1_amplitudes(:,1))
legend(hAx1, legend_label_f2, 'location' , 'northeastoutside');
xlabel('F1 stimulus amplitude [\mum]')
ylabel('logSNR at 23Hz [-]')

% Plot 200Hz responses
hAx2 = subplot(1, 2, 2);
hold on
% Use different color for different stim f2 amps.
for f1_amp_id = 1:n_f1_amp
    x_values = 1:n_f2_amp;
    %x_values = x_values + h_offset(f1_amp_id);
    errb = errorbar(x_values, mean_respf2_matrix(f1_amp_id, :), ...
             stderr_respf2_matrix(f1_amp_id, :), ...
             'color', cmp(f1_amp_id, :));
    % Remove legend for error bar
    errb_annotation = get(errb, 'Annotation');
    set(get(errb_annotation,'LegendInformation'),...
            'IconDisplayStyle','off');
         
    scatt = scatter(x_values, mean_respf2_matrix(f1_amp_id, :), ...
                    6^2, cmp(f1_amp_id, :), 'filled');
    % Remove legend for scatter
    %scatt_annotation = get(scatt, 'Annotation');
    %set(get(scatt_annotation,'LegendInformation'),...
    %        'IconDisplayStyle','off');
end
hold off
xlim([0 n_f2_amp+1])
xticks(1:n_f2_amp)
xticklabels(stim_f2_amplitudes(1,:))
legend(hAx2, legend_label_f1, 'location' , 'northeastoutside');
xlabel('F2 stimulus amplitude [\mum]')
ylabel('logSNR at 200Hz [-]')

set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [5, 5, x_width, y_width/2];
if img_fmt == "-depsc" || img_fmt == "-dpdf"
    print(gcf, img_fmt, ['FigRev1_RespAsFunc_onech']);
elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
    print(gcf, img_fmt, ['FigRev1_RespAsFunc_onech'], '-r300');
end
%}
%% Plot  all bipolar channels for each sessions.
% Take mean across trials first. 
% Then, take mean across channels
% Note that # of stimulus amplitude are different across sessions.

% Find files
% find list of names
cat_names = importdata(fullfile(data_path, 'raw', 'raw_data_filenames_abbre.txt'));
% 6x6 stimulus conditions
cat_group1 = {cat_names{1:2}};
% 4x4 stimulus conditions
cat_group2 = {cat_names{4:6}}; 
% 5x5 stimulus condtitions
cat_group3 = {cat_names{3}, cat_names{7:9}};

for cat_group_id = 1:3
    switch cat_group_id
        case 1
            cat_group = cat_group1;
        case 2
            cat_group = cat_group2;
        case 3
            cat_group = cat_group3;
    end % switch cat_group_id 

    % Plot responses as a function of stimulus amplitudes
    figure();
    % Loop through areas
    for area = 1:2    
        for cat_id = 1:length(cat_group)
            cat_name = cat_group{cat_id};
            file_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'; %evkdpwr

            % Find files
            loadnames = dir(fullfile(data_dir, cat_name, file_type, '*.mat'));
            nCond = numel(loadnames)/2; % number of conditions per area 

            % Loop through conditions in the session.
            for condid = 1:nCond

                % Load logSNR/VELogP data.
                fprintf('Loading area %i/%i, loading data %2i/%i\n', ...
                    area, 2, condid, nCond)
                load(fullfile(data_dir, cat_name, file_type, ...
                    loadnames((area-1)*nCond + condid).name))

                % if it's the first condition, 
                % 1) check whether data included unipolar, and
                % 2) preallocate frequencies of interest (per area).
                % 3) preallocate stimulus amplitude.
                if condid==1
                    % 1) check if unipolar id included in srn data. If included,
                    %   unipolar should be removed. 
                    if strcmp(data.label{1}(1:3), 'raw')
                        % the first bipolar channel
                        bipolarchid = 1 + prod(data.custom.spatialconfig);
                        % # of bipolar channels
                        nChan = data.custom.nsignals - bipolarchid + 1; 
                    else
                        bipolarchid = 1;
                        nChan = data.custom.nsignals;
                    end

                    % 2) preallocate data
                    nTrials = size(data.trial, 3);
                    resp_matrix_eacharea = zeros(nChan, length(foi), ...
                                        data.custom.subplotconfig(1),...
                                        data.custom.subplotconfig(2),...
                                        nTrials);

                    % 3) preallocate stim info matrix
                    stim_f1_amplitudes = zeros(size(data.custom.subplotconfig));
                    stim_f2_amplitudes = zeros(size(data.custom.subplotconfig));
                    % Get id and actual value of foi in continuous data
                    [foivals, foiids] = find_closest(data.freq{1}, foi);

                else % condid>1
                    % check that the value of frequency is right
                    if ~isequal(data.freq{1}(foiids(1)), foivals(1))
                        warning(['Frequency indicies inconsistent.', ... 
                            ' Recalculating. (a=%i, c=%i)'], area, condid)

                        [foivals, foiids] = find_closest(data.freq{1}, foi);
                    end
                end % if c=1

                % Set amplitude id
                f1_id = ceil(condid / data.custom.subplotconfig(2));
                f2_id = mod(condid, data.custom.subplotconfig(2));
                if f2_id == 0
                    f2_id = data.custom.subplotconfig(2);
                end % f2_ind == 0
                % get amplitude information from file name
                stim_f1_amplitudes(f1_id, f2_id) = str2double(loadnames((area-1)*nCond + condid).name(23:25));
                stim_f2_amplitudes(f1_id, f2_id) = str2double(loadnames((area-1)*nCond + condid).name(32:34));

                % data.trial = channels x frequencies x trials 
                data_bipolar = data.trial(bipolarchid:end,:,:);

                % Extract responses at foi (mean across trials)
                % channels x freqs x F1 stim amp x F2 stim amp x trials
                resp_matrix_eacharea(:, :, f1_id, f2_id , :) = data_bipolar(:, foiids, :);

            end % condid=1:nCond
            
            % mean/stderr resp at f1 and f2 across trials for each channel
            mean_respf1 = ...
                squeeze(mean(resp_matrix_eacharea(:, (foi==f1), :, :, :), 5));
            stderr_respf1 = ...
                squeeze(std(resp_matrix_eacharea(:, (foi==f1), :, :, :), 0, 5))... 
                / sqrt(nTrials);

            mean_respf2 = ...
                squeeze(mean(resp_matrix_eacharea(:, (foi==f2), :, :, :), 5));
            stderr_respf2 = ...
                squeeze(std(resp_matrix_eacharea(:, (foi==f2), :, :, :), 0, 5)) ...
                / sqrt(nTrials);
            
            % Concatenate matrix across sessions along channels
            if cat_id == 1
                mean_respf1_matrix = mean_respf1;
                stderr_respf1_matrix = stderr_respf1;
                mean_respf2_matrix = mean_respf2;
                stderr_respf2_matrix = stderr_respf2;
            else
                mean_respf1_matrix = ...
                    cat(1, mean_respf1_matrix, mean_respf1);
                stderr_respf1_matrix = ...
                    cat(1, stderr_respf1_matrix, stderr_respf1);
                mean_respf2_matrix = ...
                    cat(1, mean_respf2_matrix, mean_respf2);
                stderr_respf2_matrix = ...
                    cat(1, stderr_respf2_matrix, stderr_respf2);
            end % if cat_id == 1
            % Keep stimulus condition of Session 1-3(C20110511_R02)
            if strcmp(cat_name, 'C20110511_R02')
                stim_f2_amplitudes_sp = stim_f2_amplitudes;
            end
        end % for cat_id = 1:length(cat_names)
        
        % mean of mean/stderr resp across channels
        mean_mean_respf1_matrix = ...
            squeeze(nanmean(mean_respf1_matrix, 1)); 
        mean_stderr_respf1_matrix = ...
            squeeze(nanmean(stderr_respf1_matrix, 1));
        mean_mean_respf2_matrix = ...
            squeeze(nanmean(mean_respf2_matrix, 1));
        mean_stderr_respf2_matrix = ...
            squeeze(nanmean(stderr_respf2_matrix, 1));

        % Plot responses as a function of stimulus amplitudes
        % Legend label
        n_f1_amp = size(stim_f1_amplitudes, 1);
        n_f2_amp = size(stim_f1_amplitudes, 2);
        legend_label_f1 = {};
        legend_label_f2 = {};        
        xlabel_f2_sp = {}; % for Session 1-3
        for f_i = 1:n_f1_amp
            legend_label_f1{f_i} = ...
                [num2str(stim_f1_amplitudes(f_i, 1)), '\mum'];
            if n_f1_amp == 5        
                legend_label_f2{f_i} = ...
                    [num2str(stim_f2_amplitudes(1, f_i)), '(',...
                     num2str(stim_f2_amplitudes_sp(1, f_i)), ') ', '\mum'];
                 xlabel_f2_sp{f_i} =  ...
                     [num2str(stim_f2_amplitudes(1, f_i)), '(',...
                     num2str(stim_f2_amplitudes_sp(1, f_i)), ')'];
            else
                legend_label_f2{f_i} = ...
                    [num2str(stim_f2_amplitudes(1, f_i)), '\mum'];
            end % if n_f1_amp == 5
        end % if n_f1_amp == 5
        % Set coloar
        cmp_tmp = colormap('parula');%'gray');
        %color_step = floor((size(cmp_tmp, 1) - 40)/n_f1_amp);
        %cmp = flipud(cmp_tmp(1:color_step:256-40, :));
        color_step = floor(size(cmp_tmp-10, 1)/n_f1_amp);
        cmp = flipud(cmp_tmp(1:color_step:256-10, :));

        % plot 23Hz responses
        fig_id = area; %1 + 2*(area - 1);
        hAx1 = subplot(2, 2, fig_id);
        hold on

        % Use different color for different stim f2 amps.
        for f2_amp_id = 1:n_f2_amp
            x_values = 1:n_f1_amp;
            %x_values = x_values + h_offset(f2_amp_id);
            errb = errorbar(x_values, mean_mean_respf1_matrix(:, f2_amp_id), ...
                     mean_stderr_respf1_matrix(:, f2_amp_id), ...
                     'color', cmp(f2_amp_id, :));
            % Remove legend for error bar
            errb_annotation = get(errb, 'Annotation');
            set(get(errb_annotation,'LegendInformation'),...
                    'IconDisplayStyle','off');

            scatt = scatter(x_values, mean_mean_respf1_matrix(:, f2_amp_id), ...
                            6^2, cmp(f2_amp_id, :), 'filled');
            % Remove legend for scatter
            %scatt_annotation = get(scatt, 'Annotation');
            %set(get(scatt_annotation,'LegendInformation'),...
            %        'IconDisplayStyle','off');
        end
        hold off
        xlim([0 n_f1_amp+1])
        xticks(1:n_f1_amp)
        xticklabels(stim_f1_amplitudes(:,1))
        legend(hAx1, legend_label_f2, 'location' , 'northeastoutside');
        xlabel('F1 stimulus amplitude [\mum]')
        ylabel('logSNR at 23Hz [dB]')

        title(['S', num2str(area)]);

        % Plot 200Hz responses
        fig_id = 2 + area; %2 + 2*(area - 1);
        hAx2 = subplot(2, 2, fig_id);
        hold on
        % Use different color for different stim f2 amps.
        for f1_amp_id = 1:n_f1_amp
            % f2 stim amp varies along 2nd dim.
            start_id = 1 + n_f2_amp * (f1_amp_id - 1);
            end_id = n_f2_amp * f1_amp_id;
            x_values = 1:n_f2_amp;
            %x_values = x_values + h_offset(f1_amp_id);
            errb = errorbar(x_values, mean_mean_respf2_matrix(f1_amp_id, :), ...
                     mean_stderr_respf2_matrix(f1_amp_id, :), ...
                     'color', cmp(f1_amp_id, :));
            % Remove legend for error bar
            errb_annotation = get(errb, 'Annotation');
            set(get(errb_annotation,'LegendInformation'),...
                    'IconDisplayStyle','off');

            scatt = scatter(x_values, mean_mean_respf2_matrix(f1_amp_id, :), ...
                            6^2, cmp(f1_amp_id, :), 'filled');
            % Remove legend for scatter
            %scatt_annotation = get(scatt, 'Annotation');
            %set(get(scatt_annotation,'LegendInformation'),...
            %        'IconDisplayStyle','off');
        end
        hold off
        xlim([0 n_f2_amp+1])
        xticks(1:n_f2_amp)

        if n_f2_amp == 5
            xticklabels(xlabel_f2_sp);
        else
            xticklabels(stim_f2_amplitudes(1,:))
        end
        legend(hAx2, legend_label_f1, 'location' , 'northeastoutside');
        xlabel('F2 stimulus amplitude [\mum]')
        ylabel('logSNR at 200Hz [dB]')

        title(['S', num2str(area)]);
    end % for area = 1:2
    % Set figure properties
    %suptitle([cat_name(1:9), ' ', cat_name(11:end)]);
    set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
    set(findall(gcf,'-property','FontName'), 'FontName', ftype);
    %set(gcf,'color','w');
    set(gcf,'renderer','Painters');
    f=gcf;
    f.Units = 'centimeters';
    f.Position = [5, 5, x_width, y_width];
    if img_fmt == "-depsc" || img_fmt == "-dpdf"
        print(gcf, img_fmt, ...
             ['FigRev1_RespAsFunc_allch_numstimcond', ...
              num2str(n_f1_amp), '_eacharea']);
    elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
        print(gcf, img_fmt, ...
              ['FigRev1_RespAsFunc_allch_numstimcond', ...
               num2str(n_f1_amp), '_eacharea'], '-r300');
    end        
end %for cat_group_id = 1:3
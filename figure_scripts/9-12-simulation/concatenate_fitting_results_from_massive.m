%% Concatenate all channels data
% Concatenate fitting results from each channel into one data strucutre

func_version = 'v2';

% data directory
% only the best setting
dir_name_best_setting = ['./rect_top10_results_best_setting_', func_version, '/'];

addpath(genpath(dir_name_best_setting));

for area_id = 1:2
    % File name     
    file_name_best_setting = ['best_setting_rect_S', num2str(area_id)];
    
    full_name = fullfile([dir_name_best_setting, file_name_best_setting, ...
                         '*.mat']);
    % Get all files
    files = dir(full_name); 
    
    for channel_id = 1:size(files, 1)
        
        file_name = ['best_setting_rect_S', num2str(area_id), '_id_', ...
                    num2str(channel_id), '.mat'];        
        
        % Load the data
        load(file_name);
        
        % Store data per 
        each_channel(channel_id).logsnr = observed_logsnr;
        each_channel(channel_id).each_model = each_model_best;
        each_channel(channel_id).session = session;
        each_channel(channel_id).bipolar_ch = bipolar_ch;       
    end % channe_id = 1:size(files, 1)
    
    each_area(area_id).area =  ['S', num2str(area_id)]; % area
    each_area(area_id).freqs = freqs; % frequencies
    each_area(area_id).each_channel = each_channel;
    clear each_channel
end

save(['Results_Rect_top10_', func_version, '.mat'], 'each_area');

%% HSq
% data directory
% only the best setting
dir_name_best_setting = ['./hsq_top10_results_best_setting_', func_version, '/'];

addpath(genpath(dir_name_best_setting));

for area_id = 1:2
    % File name     
    file_name_best_setting = ['best_setting_hsq_S', num2str(area_id)];
    
    full_name = fullfile([dir_name_best_setting, file_name_best_setting, ...
                         '*.mat']);
    % Get all files
    files = dir(full_name); 
    
    for channel_id = 1:size(files, 1)
        
        file_name = ['best_setting_hsq_S', num2str(area_id), '_id_', ...
                    num2str(channel_id), '.mat'];
        
        % Load the data
        load(file_name);
        
        % Store data per 
        each_channel(channel_id).logsnr = observed_logsnr;
        each_channel(channel_id).each_model = each_model_best;
        each_channel(channel_id).session = session;
        each_channel(channel_id).bipolar_ch = bipolar_ch;       
    end % channe_id = 1:size(files, 1)
    
    each_area(area_id).area =  ['S', num2str(area_id)]; % area
    each_area(area_id).freqs = freqs; % frequencies
    each_area(area_id).each_channel = each_channel;
    clear each_channel
end

save(['Results_HSq_top10_', func_version, '.mat'], 'each_area');

%{
% Save data 
% all settings 
dir_name_all_settings = './rect_top10_results_all_settings/';

file_name_all_settings = ['all_settings_rect_S', num2str(area_id), '_id_', ...
             num2str(channel_id), '.mat'];
save([dir_name_all_settings, file_name_all_settings], 'each_model_all', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');

% only the best setting
dir_name_best_setting = './rect_top10_results_best_setting/';
save([dir_name_best_setting, file_name_best_setting], 'each_model_best', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');
%}
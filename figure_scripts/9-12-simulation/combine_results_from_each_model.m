function combine_results_from_each_model(channel_id, area_id)
%%Combine results from each model to one  
% Error = Sum across fundamentals, harmonics, and intermodulations.
% Input:    channel_id = bipolar channel id. (each area)
%           area_id = 1 or 2 (i.e. S1 or S2)
% Output:   Save fitting results.


%% Parameters
% channel_id = 99;
% area_id = 1;

% Add path
addpath(genpath('../../in_path'));
addpath(genpath('../../scripts'));

%% Load results from each model
% all settings 
dir_name_all_settings_each = './hsq_top10_results_all_settings_v2_each/';
% only the best setting
dir_name_best_settings_each = './hsq_top10_results_best_setting_v2_each/';

for model_id = 1:4
    
    % Load results from a model (all results file)
    file_name_all_settings = ...
        ['all_settings_hsq_S', num2str(area_id), '_id_', ...
         num2str(channel_id), '_model_id_', num2str(model_id), '.mat'];
    loaded_data = ...
        load([dir_name_all_settings_each, file_name_all_settings]);
    
    % Store results from all settings 
    each_model_all(model_id).error = ...
        loaded_data.each_model_all(model_id).error;
    each_model_all(model_id).parameters = ...
        loaded_data.each_model_all(model_id).model_signal;
    each_model_all(model_id).model_signal = ...
        loaded_data.each_model_all(model_id).model_signal;
    each_model_all(model_id).model_logsnr = ...
        loaded_data.each_model_all(model_id).model_logsnr;
    
    % Load best result file
    file_name_best_setting = ...
        ['best_setting_hsq_S', num2str(area_id), '_id_', ...
         num2str(channel_id), '_model_id_', num2str(model_id), '.mat'];
    loaded_best_data = ...
        load([dir_name_best_settings_each, file_name_best_setting]);
    % Store the best across all settings
    each_model_best(model_id).error = ...
        loaded_best_data.each_model_best(model_id).error;
    each_model_best(model_id).parameters = ...
        loaded_best_data.each_model_best(model_id).parameters;
    each_model_best(model_id).model_signal = ...
        loaded_best_data.each_model_best(model_id).model_signal;
    each_model_best(model_id).model_logsnr = ...
        loaded_best_data.each_model_best(model_id).model_logsnr;   
end

% Rename variable for storing them
observed_logsnr = loaded_data.observed_logsnr;
freqs = loaded_data.freqs;
session = loaded_data.session;
bipolar_ch = loaded_data.bipolar_ch;
area_id = loaded_data.area_id;

% Save data 
% all settings 
dir_name_all_settings = './hsq_top10_results_all_settings_v2/';

if ~exist(dir_name_all_settings, 'dir')
    mkdir(dir_name_all_settings);
end % ~exist(dir_name_all_settings, 'dir');

file_name_all_settings = ['all_settings_hsq_S', num2str(area_id), '_id_', ...
             num2str(channel_id), '.mat'];
save([dir_name_all_settings, file_name_all_settings], 'each_model_all', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');         

% only the best setting
dir_name_best_setting = './hsq_top10_results_best_setting_v2/';
if ~exist(dir_name_best_setting, 'dir')
    mkdir(dir_name_best_setting);
end % ~exist(dir_name_all_settings, 'dir');
file_name_best_setting = ['best_setting_hsq_S', num2str(area_id), '_id_', ...
             num2str(channel_id), '.mat'];
save([dir_name_best_setting, file_name_best_setting], 'each_model_best', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');

end
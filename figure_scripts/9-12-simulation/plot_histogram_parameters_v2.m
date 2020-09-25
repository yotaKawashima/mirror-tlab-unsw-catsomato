%% Histogram of parameters
% Plot histograms of best parameters for each model across channels
% model 
% Show only 1%~99% data

%% Rect compare S1 and S2 for each model
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
Rect_data = load('Results_Rect_top10_v2.mat');

% Get parameters from each area
num_params = [2, 3, 3, 4];
name_params = ['a', 'b', 'c', 'd'];
disp('Rect');
subplot_ids = [1, 2, 5, 6, 7, 9, 10, 11, 13, 14, 16, 15];

for area_id = 1:2
    subplot_incre = 1;
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for model_id = [1, 3, 2, 4]
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id = 1:size(best_parameters, 2)
            subplot(4, max(num_params), subplot_ids(subplot_incre));
            % Get 5%~95%
            best_parameters_now = best_parameters(:, param_id);
            trim_threshold = prctile(best_parameters_now,  [3,97]);
            % Get index
            lower_ = best_parameters_now>trim_threshold(1);
            upper_ = best_parameters_now<trim_threshold(2);
            % Remove data outside of the range
            trimed_parameters = best_parameters_now(logical(lower_ .* upper_));
            plot_line = cdfplot(trimed_parameters);
            set(plot_line, 'linewidth', 3);
            xlabel('parameter [-]');
            ylabel('probability [-]');
            title(['S', num2str(area_id), ', model id ', num2str(model_id),...
                   ', parameter ', name_params(param_id)]);
           
           subplot_incre = subplot_incre + 1;
    
        end % for pram_id         
    end % for model_id 
end %for area_id


%% Sq 
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
Sq_data = load('Results_Sq_top10_v2.mat');

% Get parameters from each area
num_params = [2, 3];
name_params = ['a', 'b', 'c'];
disp('Sq');
subplot_ids = [1, 2, 4, 5, 6];

for area_id = 1:2
    subplot_incre = 1;
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for model_id = [1, 2]
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id = 1:size(best_parameters, 2)
            subplot(2, max(num_params), subplot_ids(subplot_incre));
            % Get 5%~95%
            best_parameters_now = best_parameters(:, param_id);
            %trim_threshold = prctile(best_parameters_now, [30 70]);
            % Get index
            lower_ = best_parameters_now>trim_threshold(1);
            upper_ = best_parameters_now<trim_threshold(2);
            % Remove data outside of the range
            trimed_parameters = best_parameters_now;%best_parameters_now(logical(lower_ .* upper_));
            plot_line = cdfplot(trimed_parameters);
            set(plot_line, 'linewidth', 3);
            xlabel('parameter [-]');
            ylabel('probability [-]');
            title(['S', num2str(area_id), ', model id ', num2str(model_id),...
                   ', parameter ', name_params(param_id)]);
            
            subplot_incre = subplot_incre + 1;
           
        end % for pram_id         
    end % for model_id 
end %for area_id

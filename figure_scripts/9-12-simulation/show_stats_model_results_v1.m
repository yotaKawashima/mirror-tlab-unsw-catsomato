%% Stats model results
% 1. Wilcoxon signed rank test to compare four rect models for each
%    somatosensory area. 
% 2. Get some percentiles of the min difference for rect models and sq
%    models.
% model 
% first part
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
% Second part
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)

Rect_data = load('Results_Rect_top10.mat');

%% Rect compare S1 and S2 for each model
% Get error from each area
error_each_area = struct();
disp('Rect');
for area_id = 1:2
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
    
    % Compute difference between errors
    difference_errors = nan(size(errors, 1), 6); % 4C2=6

    % model1 - model3:{Rect(X)+Rect(Y)} - {Rect(X)+Rect(Y)+Rect(XY)}
    difference_errors(:, 1) = errors(:, 1) - errors(:, 3);    
    % model1 - model2:{Rect(X)+Rect(Y)} - {Rect(X)+Rect(Y)+Rect(X)Rect(Y)}
    difference_errors(:, 2) = errors(:, 1) - errors(:, 2);    
    % model1 - model4:{Rect(X)+Rect(Y)} - {Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)}
    difference_errors(:, 3) = errors(:, 1) - errors(:, 4);    
    % model3 - model2:{Rect(X)+Rect(Y)+Rect(XY)} - {Rect(X)+Rect(Y)+Rect(X)Rect(Y)}
    difference_errors(:, 4) = errors(:, 3) - errors(:, 2);        
    % model3 - model4:{Rect(X)+Rect(Y)+Rect(XY)} - {Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)}
    difference_errors(:, 5) = errors(:, 3) - errors(:, 4);   
    % model2 - model4:{Rect(X)+Rect(Y)+Rect(X)Rect(Y)} - {Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)}
    difference_errors(:, 6) = errors(:, 2) - errors(:, 4);        
     
    % Store error data for each area
    error_each_area(area_id).errors = errors;    
    error_each_area(area_id).difference_errors = difference_errors;  
    
    % Show some percentiles
    disp('percentiles: ');        
    for error_id = [1, 3, 2, 4]
        disp(['model' num2str(error_id)]);
        disp(['min ', num2str(min(errors(:, error_id)))]);
        disp(['max ', num2str(max(errors(:, error_id)))]);
        disp(['25th ', num2str(prctile(errors(:, error_id), 25))]);
        disp(['50th ', num2str(prctile(errors(:, error_id), 50))]);
        disp(['75th ', num2str(prctile(errors(:, error_id), 75))]);        
    end
end

% Wilcoxon singed rank test for each somatosensory area.
% Plot diff distribution for each area
rect_ttest_stats = struct();

for area_id = 1:2
    difference_errors = error_each_area(area_id).difference_errors;
    
    % Across all comparison
    h_wilcoxon_test = [];
    p_wilcoxon_test = [];
    stats_wilcoxon_test = [];
    for comparison_id = 1:6
        
        % wilcoxon signed rank test  
        % Whether greater than 0 or not
        [p_w, h_w, stats_w] = signrank(difference_errors(:, comparison_id),...
                                          0, 'Tail', 'right'); 
        
        p_wilcoxon_test = [p_wilcoxon_test; p_w];
        h_wilcoxon_test = [h_wilcoxon_test; h_w];
        stats_wilcoxon_test = [stats_wilcoxon_test; stats_w];
    end
    
    % Show stats results 
    zvals = [stats_wilcoxon_test.zval]';
    signedranks = [stats_wilcoxon_test.signedrank]';
    wilcoxon_test_stats = table(h_wilcoxon_test, p_wilcoxon_test, zvals, signedranks);
    disp(['Rect: S', num2str(area_id)]);
    disp(wilcoxon_test_stats);
    
    rect_ttest_stats(area_id).p = p_wilcoxon_test;
    rect_ttest_stats(area_id).h = h_wilcoxon_test;
    rect_ttest_stats(area_id).stats = stats_wilcoxon_test;
            
end %area_id = 1:2

%% Sq 
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
Sq_data = load('Results_Sq_top10.mat');

% Get error from each area
error_each_area = struct();
disp('Sq');
for area_id = 1:2
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);    
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        % Get the first and the third results
        % 1st: Sq(X)+Sq(Y)
        % 2nd: Sq(X)+Sq(Y)+Sq(X)Sq(Y)
        each_error_1_3 = [each_error(1), each_error(3)];
        errors = [errors; each_error_1_3];
    end % for ch_id = 1:num_channel

    % Store error data for each area
    error_each_area(area_id).errors = errors;    
    
    
    % Show some percentiles
    disp('percentiles: ');        
    for error_id = [1, 2]
        disp(['model' num2str(error_id)]);
        disp(['min ', num2str(min(errors(:, error_id)))]);
        disp(['max ', num2str(max(errors(:, error_id)))]);
        disp(['25th ', num2str(prctile(errors(:, error_id), 25))]);
        disp(['50th ', num2str(prctile(errors(:, error_id), 50))]);
        disp(['75th ', num2str(prctile(errors(:, error_id), 75))]);
    end

end %area_id = 1:2


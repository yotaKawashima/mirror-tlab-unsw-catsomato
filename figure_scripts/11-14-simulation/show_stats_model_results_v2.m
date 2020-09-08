%% Stats model results
% 1. Get some percentiles of the min difference for rect models and sq
%    models.

%% Rect
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)

Rect_data = load('Results_Rect_top10_v2.mat');
% Get error from each area
error_each_area = struct();
num_params = [2, 3, 3, 4];
disp('Rect');
for area_id = 1:2
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    disp('------');

    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
     
    % Store error data for each area
    error_each_area(area_id).errors = errors;    
    
    AICs = nans(1,4);
    % Show some percentiles
    disp('percentiles: ');        
    model_ids = [1, 3, 2, 4];
    for order_id = 1:4
        model_id = model_ids(order_id);
        disp(['---model' num2str(model_id)]);
        disp(['min ', num2str(min(errors(:, model_id)))]);
        disp(['max ', num2str(max(errors(:, model_id)))]);
        disp(['25th ', num2str(prctile(errors(:, model_id), 25))]);
        disp(['50th ', num2str(prctile(errors(:, model_id), 50))]);
        disp(['75th ', num2str(prctile(errors(:, model_id), 75))]);       
        
        % Compute AIC
        sigma = sum(errors(:, model_id).^2) / num_channels;
        k_val = num_params(model_id)+ 1;
        AIC = num_channels*log(sigma) + 2*k_val; % base = e 
        AICs(1, order_id) = AIC;
        disp(['sigma ', num2str(sigma)]);
        disp(['k ', num2str(k_val)]);
        disp(['AIC ', num2str(AIC)]);
    end % for model_id =
    disp('------');

end


%% Sq 
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
Sq_data = load('Results_Sq_top10_v2.mat');

% Get error from each area
error_each_area = struct();
num_params = [2, 3];
disp('Sq');
for area_id = 1:2
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);  
    disp('------');
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        % Get the first and the third results
        % 1st: Sq(X)+Sq(Y)
        % 2nd: Sq(X)+Sq(Y)+Sq(X)Sq(Y)
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel

    % Store error data for each area
    error_each_area(area_id).errors = errors;    
    
    
    % Show some percentiles
    disp('percentiles: ');        
    for model_id = [1, 2]
        disp(['---model' num2str(model_id)]);
        disp(['min ', num2str(min(errors(:, model_id)))]);
        disp(['max ', num2str(max(errors(:, model_id)))]);
        disp(['25th ', num2str(prctile(errors(:, model_id), 25))]);
        disp(['50th ', num2str(prctile(errors(:, model_id), 50))]);
        disp(['75th ', num2str(prctile(errors(:, model_id), 75))]);

        % Compute AIC
        sigma = sum(errors(:, model_id).^2) / num_channels;
        k_val = num_params(model_id)+ 1;
        AIC = num_channels*log(sigma) + 2*k_val; % base = e 
        disp(['sigma ', num2str(sigma)]);
        disp(['k ', num2str(k_val)]);
        disp(['AIC ', num2str(AIC)]);
    end % model_id = 
    disp('------');
end %area_id = 1:2


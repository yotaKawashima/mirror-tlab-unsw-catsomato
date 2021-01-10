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
    
    AICs = nan(1,4);
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
    
    % Rect
    rect_aic(area_id).aic = AICs;

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
    
    AICs = nan(1,2);
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
        AIC = num_channels*log(sigma) + 2*k_val; % base = e 
        AICs(1, model_id) = AIC;
        disp(['sigma ', num2str(sigma)]);
        disp(['k ', num2str(k_val)]);
        disp(['AIC ', num2str(AIC)]);
    end % model_id = 
    
    % Rect
    sq_aic(area_id).aic = AICs;

    disp('------');
end %area_id = 1:2


%% Plot AIC
colours = [[0, 0.4470, 0.7410];...
           [0.8500, 0.3250, 0.0980]];
       
shift_value = 7;       
figure(); clf
hold on
xs = [6,5,4,3,2,1] + shift_value;
    
p = gobjects(1,2);
for area_id = 1:2
    rect_aic_data = rect_aic(area_id).aic;
    sq_aic_data = sq_aic(area_id).aic;
    
    concatenated_aic = [rect_aic_data, sq_aic_data];
    if area_id == 2
        xs = xs - shift_value;
    end
    barh(xs, concatenated_aic, 'FaceColor', colours(area_id, :));
    %p(area_id) = plot(concatenated_aic, xs, 'Color', 'k' , 'LineWidth', 3);
    
end
hold off

xlim([500, 1500]);
% Remove legend for lines
%p_annotation = get(p, 'Annotation');
%for i_pline = 1:length(p)
%    set(get(p_annotation{i_pline,1},'LegendInformation'),...
%        'IconDisplayStyle','off');
%end % i_vline = 1:length(v)

xlabel('AIC');
y_labels = {'Rect(X)+Rect(Y)';'Rect(X)+Rect(Y)+Rect(XY)';...
           'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';...
           'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)';...
           'Sq(X)+Sq(Y)';'Sq(X)+Sq(Y)+Sq(XY)';...
           'Rect(X)+Rect(Y)';'Rect(X)+Rect(Y)+Rect(XY)';...
           'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';...
           'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)';...
           'Sq(X)+Sq(Y)';'Sq(X)+Sq(Y)+Sq(XY)'};
    
yticks(flip([xs+shift_value, xs]));
yticklabels(flip(y_labels));
legend({'S1', 'S2'});

set(gca, 'fontsize', 16);
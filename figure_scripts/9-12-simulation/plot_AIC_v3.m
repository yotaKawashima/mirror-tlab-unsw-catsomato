%% Stats model results
% 1. Get some percentiles of the min difference for rect models and sq
%    models.

%% Rect
% aRect(X) + bRect(Y)                                   model_id = 1
% aRect(X) + bRect(Y) + cRect(X)Rect(Y)                 model_id = 2
% aRect(X) + bRect(Y) + cSRect(XY)                      model_id = 3
% aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY)     model_id = 4

img_fmt = '-dpng';

fsize = 8.5; % font size
ftype = 'Arial'; % font type
x_width = 13; % fig width
y_width = 5; % fig height


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

%% HSq
% aHSq(X) + bHSq(Y)                                model_id = 1
% aHSq(X) + bHSq(Y) + cHSq(X)HSq(Y)                model_id = 2
% aHSq(X) + bHSq(Y) + cHSq(XY)                     model_id = 3
% aHSq(X) + bHSq(Y) + cHSq(X)HSq(Y) + dHSq(XY)     model_id = 4
HSq_data = load('Results_HSq_top10_v2.mat');

% Get error from each area
error_each_area = struct();
num_params = [2, 3, 3, 4];
disp('HSq');
for area_id = 1:2
    data_each_channel = HSq_data.each_area(area_id).each_channel;
    
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
    
    % HSq
    hsq_aic(area_id).aic = AICs;

    disp('------');

end %for area_id = 1:2
%% Plot AIC
colours = [[0, 0.4470, 0.7410];...
           [0.8500, 0.3250, 0.0980]];

xs = [8,7,6,5,4,3,2,1];
    
for area_id = 1:2
    figure(); clf

    rect_aic_data = rect_aic(area_id).aic;
    hsq_aic_data = hsq_aic(area_id).aic;
    
    concatenated_aic = [rect_aic_data, hsq_aic_data];
    barh(xs, concatenated_aic, 'FaceColor', colours(area_id, :));
    %p(area_id) = plot(concatenated_aic, xs, 'Color', 'k' , 'LineWidth', 3);
    
    y_labels = {'Rect(X)+Rect(Y)';...
                'Rect(X)+Rect(Y)+Rect(XY)';...
                'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';...
                'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)';...
                'HSq(X)+HSq(Y)';...
                'HSq(X)+HSq(Y)+HSq(XY)';...
                'HSq(X)+HSq(Y)+HSq(X)HSq(Y)';...
                'HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)';...
               };
    xlabel('AIC');
    yticks(flip(xs));
    yticklabels(flip(y_labels));
    if area_id == 1
        xlim([1200, 1350]);
        title('S1');
    elseif area_id == 2
        xlim([750, 900]);      
        title('S2');        
    end
    
    set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
    set(findall(gcf,'-property','FontName'), 'FontName', ftype);
    set(gcf,'renderer','Painters');
    f=gcf;
    f.Units = 'centimeters';
    f.Position = [10, 10, x_width, y_width];
    % Print
    switch area_id 
        case 1
            filename = 'figure12_a';
        case 2             
            filename = 'figure12_b';
    end
    if img_fmt == "-depsc" || img_fmt == "-dpdf"   
        print(gcf, img_fmt, filename);
    elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
        print(gcf, img_fmt, filename, '-r300');
    end

end



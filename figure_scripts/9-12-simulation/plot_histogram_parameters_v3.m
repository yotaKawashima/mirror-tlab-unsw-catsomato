%% Histogram of parameters
% Plot histograms of best parameters for each model across channels
% model 


%% Rect compare S1 and S2 for each model
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
Rect_data = load('Results_Rect_top10_v2.mat');

img_fmt = "-dpdf";

fsize = 8.5; % font size
ftype = 'Arial'; % font type
x_width = 13; % fig width
y_width = 5; % fig height


% Get parameters from each area
num_params = [2, 3, 3, 4];
name_params = ['a', 'b', 'c', 'd'];
disp('Rect');

model_orders = [1, 3, 2, 4];

%blue, red, puple, yellow
%simplest, Rect(XY), Rect(X)Rect(Y), most complex
colours = [[0, 0.4470, 0.7410];...
           [0.8500, 0.3250, 0.0980];...
           [0.4940, 0.1840, 0.5560];...
           [0.9290, 0.6940, 0.1250]];

for area_id = 1:2
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for order_id = 1:4
        model_id = model_orders(order_id);
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id = 1:size(best_parameters, 2)
            subplot(2, 2, param_id);
            hold on
            % Get 5%~95%
            best_parameters_now = best_parameters(:, param_id);
            trim_threshold = prctile(best_parameters_now,  [3,97]);
            % Get index
            lower_ = best_parameters_now>trim_threshold(1);
            upper_ = best_parameters_now<trim_threshold(2);
            % Remove data outside of the range
            trimed_parameters = best_parameters_now(logical(lower_ .* upper_));
            plot_line = cdfplot(trimed_parameters);
            set(plot_line, 'linewidth', 1);
            set(plot_line, 'color', colours(order_id, :));
            xlabel('parameter [-]');
            ylabel('probability [-]');
            title(['S', num2str(area_id), ', parameter ', ...
                   name_params(param_id)]);
            hold off    
            if param_id == 1 && model_id == 4
                legend({'Rect(X)+Rect(Y)', 'Rect(X)+Rect(Y)+Rect(XY)',...
                        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)', ...
                        'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)'});
            end % if param_id == 1
           
        end % for pram_id   
    end % for model_id   
    
    set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
    set(findall(gcf,'-property','FontName'), 'FontName', ftype);
    set(gcf,'renderer','Painters');
    f=gcf;
    f.Units = 'centimeters';
    f.Position = [10, 10, x_width, y_width*2];
    filename = ['Sup_figure6_S', num2str(area_id)];
    % Print
    if img_fmt == "-depsc" || img_fmt == "-dpdf"   
        print(gcf, img_fmt, filename);
    elseif img_fmt == "-dtiff"
        print(gcf, img_fmt, filename, '-r300');
    end

end %for area_id



%% Sq 
%{
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
Sq_data = load('Results_Sq_top10_v2.mat');

% Get parameters from each area
num_params = [2, 3];
name_params = ['a', 'b', 'c'];
disp('Sq');

colours = [[0, 0.4470, 0.7410];...
           [0.8500, 0.3250, 0.0980]];

for area_id = 1:2
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for model_id = 1:2
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id = 1:size(best_parameters, 2)
            subplot(1, 3, param_id);
            hold on
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
            set(plot_line, 'color', colours(model_id, :));
            xlabel('parameter [-]');
            ylabel('probability [-]');
            title(['S', num2str(area_id), ', parameter ', ...
                   name_params(param_id)]);
            hold off    
            
            if param_id == 1 && model_id == 2
                legend({'model 1', 'model 2'});
            end % if param_id == 1

        end % for pram_id   
    end % for model_id     
end %for area_id
%}
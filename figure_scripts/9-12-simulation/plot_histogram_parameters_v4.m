%% Histogram of parameters
% Plot histograms of best parameters for each model across channels
% model 

%% Rect compare S1 and S2 for each model
% aRect(X) + bRect(Y)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y)
% aRect(X) + bRect(Y) + cRect(XY)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY)
Rect_data = load('Results_Rect_top10_v2.mat');

img_fmt = "-dpng";

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
        for param_id_f = 1:size(best_parameters, 2)
            % all model: aF(X)+bF(Y)+cF(X)F(Y)+dF(XY)
            % Flip c and d
            if model_id == 4
                if param_id_f == 3
                    param_id = 4;
                elseif param_id_f == 4
                    param_id = 3;
                else
                    param_id = param_id_f;
                end % if param_id_f
            else
                param_id = param_id_f;
            end % if model_id
            subplot(2, 2, param_id_f);
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
                   name_params(param_id_f)]);
            hold off    
            if param_id_f == 1 && model_id == 4
                legend({'Rect(X)+Rect(Y)', 'Rect(X)+Rect(Y)+Rect(XY)',...
                        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)', ...
                        'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)'});
            end % if param_id == 1
           
        end % for pram_id_f   
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
    elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
        print(gcf, img_fmt, filename, '-r300');
    end

end %for area_id



%% HSq 
% aHSq(X) + bHSq(Y)
% aHSq(X) + bHSq(Y) + cHSq(X)HSq(Y)
%...
HSq_data = load('Results_HSq_top10_v2.mat');
% Get parameters from each area
num_params = [2, 3, 3, 4];
name_params = ['a', 'b', 'c', 'd'];
disp('HSq');

model_orders = [1, 3, 2, 4];

%blue, red, puple, yellow
%simplest, Rect(XY), Rect(X)Rect(Y), most complex
colours = [[0, 0.4470, 0.7410];...
           [0.8500, 0.3250, 0.0980];...
           [0.4940, 0.1840, 0.5560];...
           [0.9290, 0.6940, 0.1250]];

for area_id = 1:2
    data_each_channel = HSq_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for order_id = 1:4
        model_id = model_orders(order_id);
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = ...
                data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id_f = 1:size(best_parameters, 2)
            % all model: aF(X)+bF(Y)+cF(X)F(Y)+dF(XY)
            % Flip c and d
            if model_id == 4
                if param_id_f == 3
                    param_id = 4;
                elseif param_id_f == 4
                    param_id = 3;
                else
                    param_id = param_id_f;
                end % if param_id_f
            else
                param_id = param_id_f;
            end % if model_id
            subplot(2, 2, param_id_f);
            hold on
            % Get 5%~95%
            best_parameters_now = best_parameters(:, param_id);
            %trim_threshold = prctile(best_parameters_now,  [3,97]);
            % Get index
            lower_ = best_parameters_now;%>trim_threshold(1);
            upper_ = best_parameters_now;%<trim_threshold(2);
            % Remove data outside of the range
            try
                trimed_parameters = best_parameters_now(logical(lower_ .* upper_));
                plot_line = cdfplot(trimed_parameters);
                set(plot_line, 'linewidth', 1);
                set(plot_line, 'color', colours(order_id, :));
                xlabel('parameter [-]');
                ylabel('probability [-]');
                title(['S', num2str(area_id), ', parameter ', ...
                       name_params(param_id_f)]);
                hold off    
                if param_id_f == 1 && model_id == 4
                    %legend({'HSq(X)+HSq(Y)', 'HSq(X)+HSq(Y)+HSq(XY)',...
                    %        'HSq(X)+HSq(Y)+HSq(X)HSq(Y)', ...
                    %        'HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)'});
                end % if param_id == 1
            catch 
                %
                print('ignore');
            end %try
        end % for pram_id_f   
    end % for model_id   
    
    set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
    set(findall(gcf,'-property','FontName'), 'FontName', ftype);
    set(gcf,'renderer','Painters');
    f=gcf;
    f.Units = 'centimeters';
    f.Position = [10, 10, x_width, y_width*2];
    filename = ['Sup_figure6_HSq_S', num2str(area_id)];
    % Print
    if img_fmt == "-depsc" || img_fmt == "-dpdf"   
        print(gcf, img_fmt, filename);
    elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
        print(gcf, img_fmt, filename, '-r300');
    end

end %for area_id

% narrow x axis
for area_id = 1:2
    data_each_channel = HSq_data.each_area(area_id).each_channel;
    num_channels = size(data_each_channel, 2);
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    
    figure(); clf
    for order_id = 1:4
        model_id = model_orders(order_id);
        best_parameters = nan(num_channels, num_params(model_id));
        for ch_id = 1:num_channels
            best_parameters(ch_id, :) = ...
                data_each_channel(ch_id).each_model(model_id).parameters;
        end % for ch_id
        
        % Plot parameters
        for param_id_f = 1:size(best_parameters, 2)
            % all model: aF(X)+bF(Y)+cF(X)F(Y)+dF(XY)
            % Flip c and d
            if model_id == 4
                if param_id_f == 3
                    param_id = 4;
                elseif param_id_f == 4
                    param_id = 3;
                else
                    param_id = param_id_f;
                end % if param_id_f
            else
                param_id = param_id_f;
            end % if model_id
            subplot(2, 2, param_id_f);
            hold on
            % Get 5%~95%
            best_parameters_now = best_parameters(:, param_id);
            trim_threshold = prctile(best_parameters_now,  [3,97]);
            % Get index
            lower_ = best_parameters_now>trim_threshold(1);
            upper_ = best_parameters_now<trim_threshold(2);
            % Remove data outside of the range
            try
                trimed_parameters = best_parameters_now(logical(lower_ .* upper_));
                plot_line = cdfplot(trimed_parameters);
                set(plot_line, 'linewidth', 1);
                set(plot_line, 'color', colours(order_id, :));
                xlabel('parameter [-]');
                ylabel('probability [-]');
                title(['S', num2str(area_id), ', parameter ', ...
                       name_params(param_id_f)]);
                hold off    
                if param_id_f == 1 && model_id == 4
                    %legend({'HSq(X)+HSq(Y)', 'HSq(X)+HSq(Y)+HSq(XY)',...
                    %        'HSq(X)+HSq(Y)+HSq(X)HSq(Y)', ...
                    %        'HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)'});
                end % if param_id == 1
            catch 
                %
                print('ignore');
            end %try
        end % for pram_id_f   
    end % for model_id   
    
    set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
    set(findall(gcf,'-property','FontName'), 'FontName', ftype);
    set(gcf,'renderer','Painters');
    f=gcf;
    f.Units = 'centimeters';
    f.Position = [10, 10, x_width, y_width*2];
    filename = ['Sup_figure6_HSq_narrow_S', num2str(area_id)];
    % Print
    if img_fmt == "-depsc" || img_fmt == "-dpdf"   
        print(gcf, img_fmt, filename);
    elseif img_fmt == "-dtiff" || img_fmt == "-dpng"
        print(gcf, img_fmt, filename, '-r300');
    end

end %for area_id
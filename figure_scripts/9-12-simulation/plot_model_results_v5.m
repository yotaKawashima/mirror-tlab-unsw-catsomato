%% Plot model results
% model 
% firsst part
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
% Second part
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)

%img_fmt = '-depsc2';
%img_fmt = '-dpng';
%img_fmt = '-dtiff';
img_fmt = '-dpdf';

fsize = 8.5; % font size
ftype = 'Arial'; % font type
x_width = 13; % fig width
y_width = 5; % fig height

bonferroni_correction = 6; % number of models we tested.

Rect_data = load('Results_Rect_top10_v2.mat');

%% Rect compare S1 and S2 for each model
% Get error from each area
model_orders = [1, 3, 2, 4];
error_each_area = struct();
disp('Rect');
for area_id = 1:2
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
    
    % Store error data for each area
    error_each_area(area_id).errors = errors;    
end

% Two-sample Kolmogorov-Smirnov test for each model between S1 and S2
h = [];
p = [];
ks2stat = [];
x_ks = [];
cdf1_ks = [];
cdf2_ks = [];

for model_id = 1:4
    %[h_, p_, ks2stat_] = kstest2(error_each_area(1).errors(:,model_id), ...
    %                             error_each_area(2).errors(:,model_id));
    [h_, p_, ks2stat_, x_, y1_, y2_] = kstest2_y(...
                                        error_each_area(1).errors(:,model_id), ...
                                        error_each_area(2).errors(:,model_id));
    
    h = [h; h_];
    p = [p; p_];
    ks2stat = [ks2stat; ks2stat_]; 
    x_ks = [x_ks; x_];
    cdf1_ks = [cdf1_ks; y1_];
    cdf2_ks = [cdf2_ks; y2_];
    
end % for model_id =1:4
ksstats_table = table(h, p, ks2stat);
disp(ksstats_table);

% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    errors = error_each_area(area_id).errors;
    if area_id == 1
        line_col = [0, 0.4470, 0.7410];
    elseif area_id == 2
        line_col = [0.8500, 0.3250, 0.0980];
    end
    
    for order_id = 1:4
        model_id = model_orders(order_id);
        error_now = errors(:,model_id);
        % Remove the worse 1% channels for a display purpose for
        % Rect(X)+Rect(Y) and Rect(X)+Rect(Y)+Rect(XY)
        if (area_id == 1 && model_id == 1) || ...
           (area_id == 1 && model_id == 3)
            error_threshold = prctile(error_now, 99);
            error_temp = error_now;
            error_now = error_temp(error_temp < error_threshold);
        end % if area_id == 1        
        hold on
        subplot(2, 2, order_id);
        % Plot cdf
        [h, stats] = cdfplot(error_now);
        set(h, 'linewidth', 1);
        set(h, 'color', line_col);
        if area_id ==2 % ADd line only after plotting S2 results
        % Add line showing ks stats (max distance)
        v = plot([x_ks(model_id), x_ks(model_id)],...
                 [cdf1_ks(model_id), cdf2_ks(model_id)],...
                 'color', 'k', 'linewidth', 1);
        % Ignore this line for legend 
        v_annotation = get(v, 'Annotation');
        set(get(v_annotation,'LegendInformation'),...
                    'IconDisplayStyle','off');
        end 
        xlabel('minimum differences');
        ylabel('probability');
        hold off
        switch model_id
            case 1
                title_now = 'Rect(X)+Rect(Y)';
            case 2
                title_now = 'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';
            case 3
                title_now = 'Rect(X)+Rect(Y)+Rect(XY)';
            case 4
                title_now = 'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)';
        end % swithch model_id   
       
        %title_now_3 = ['p = ', ...
        %    num2str(ksstats_table.p(model_id)*bonferroni_correction)];

        %title({title_now_1, title_now_2, title_now_3});
        title(title_now);
        
        xlim([40, 120]);
        
    end % for model_id = 1:4
end % area_id = 1:2

legend({'S1', 'S2'});
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width*2];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'figure11_a');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'figure11_a', '-r300');
end

%% Rect compare models for each area
% Plot culmtive probability distribution
colours = [[0, 0.4470, 0.7410];...
           [0.4940, 0.1840, 0.5560];...
           [0.8500, 0.3250, 0.0980];...
           [0.9290, 0.6940, 0.1250]];

figure(); clf
for area_id = 1:2
    % Get error data
    errors = error_each_area(area_id).errors;
    
    subplot(1,2, area_id);
    hold on 
    for order_id = 1:4
        model_id = model_orders(order_id); 
        error_now = errors(:, model_id);   
        % Remove the worse 1% channels for a display purpose for
        % Rect(X)+Rect(Y) and Rect(X)+Rect(Y)+Rect(XY)
        if (area_id == 1 && model_id == 1) || ...
           (area_id == 1 && model_id == 3)
            error_threshold = prctile(error_now, 99);
            error_temp = error_now;
            error_now = error_temp(error_temp < error_threshold);
        end % if area_id == 1        
        
        [h, stats] = cdfplot(error_now);
        set(h, 'linewidth', 1);
        set(h, 'color', colours(model_id,:));         
        if area_id == 1
            title('S1');
        elseif area_id == 2
            title('S2');
        end % if area_id ==1
    end % for model_id = 1:4
    xlim([40, 120]);
    xlabel('minimum differences');
    ylabel('probability');
    hold off
    
end %area_id = 1:2
legend({'Rect(X)+Rect(Y)', 'Rect(X)+Rect(Y)+Rect(XY)', ...
        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)', ...
        'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)'});
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'figure11_b');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'figure11_b', '-r300');
end

%% Sq 
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(XY)
Sq_data = load('Results_Sq_top10_v2.mat');

%% Sq compare S1 and S2 for each model
% Get error from each area
error_each_area = struct();
disp('Sq');
for area_id = 1:2
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    disp(['area S',  num2str(area_id)]);
    disp(['# of ch:', num2str(num_channels)]);    
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel

    % Store error data for each area
    error_each_area(area_id).errors = errors;    
    
end %area_id = 1:2

% Two-sample Kolmogorov-Smirnov test for each model between S1 and S2
h = [];
p = [];
ks2stat = [];
x_ks = [];
cdf1_ks = [];
cdf2_ks = [];

for model_id = 1:2
    %[h_, p_, ks2stat_] = kstest2(error_each_area(1).errors(:,model_id), ...
    %                             error_each_area(2).errors(:,model_id));
    [h_, p_, ks2stat_, x_, y1_, y2_] = kstest2_y(...
                                        error_each_area(1).errors(:,model_id), ...
                                        error_each_area(2).errors(:,model_id));
    
    h = [h; h_];
    p = [p; p_];
    ks2stat = [ks2stat; ks2stat_]; 
    x_ks = [x_ks; x_];
    cdf1_ks = [cdf1_ks; y1_];
    cdf2_ks = [cdf2_ks; y2_];
    
end % for model_id =1:4
ksstats_table = table(h, p, ks2stat);
disp(ksstats_table);

% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    errors = error_each_area(area_id).errors;
    if area_id == 1
        line_col = [0, 0.4470, 0.7410];
    elseif area_id == 2
        line_col = [0.8500, 0.3250, 0.0980];
    end
    
    for model_id = 1:2
        hold on
        subplot(1, 2, model_id);
        error_now = errors(:, model_id);
        % Plot cdf
        [h, stats] = cdfplot(error_now);
        set(h, 'linewidth', 1);
        set(h, 'color', line_col);
        if area_id ==2 % ADd line only after plotting S2 results
        % Add line showing ks stats (max distance)
        v = plot([x_ks(model_id), x_ks(model_id)],...
                 [cdf1_ks(model_id), cdf2_ks(model_id)],...
                 'color', 'k', 'linewidth', 1);
        % Ignore this line for legend 
        v_annotation = get(v, 'Annotation');
        set(get(v_annotation,'LegendInformation'),...
                    'IconDisplayStyle','off');
        end
        
        xlabel('minimum differences');
        ylabel('probability');
        hold off
        
        switch model_id
            case 1
                title_now = 'Sq(X)+Sq(Y)'; 
            case 2
                title_now = 'Sq(X)+Sq(Y)+Sq(XY)';
        end % swithch model_id
        %title_now_3 = ['p = ', ...
        %    num2str(ksstats_table.p(model_id)*bonferroni_correction)];

        %title({title_now_1, title_now_2, title_now_3});
        title(title_now);
        
        xlim([60, 200]);
        xticks([60, 100, 140, 180]);
    
    end % for model_id = 1:4

end % area_id = 1:2
legend({'S1', 'S2'});
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'Sup_figure5_ap');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'Sup_figure5_ap', '-r300');
end

%% Sq compare models for each area
% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    % Get error data
    errors = error_each_area(area_id).errors;
    
    subplot(1,2, area_id);
    hold on 
    
    for model_id = 1:2
        error_now = errors(:, model_id);        
        [h, stats] = cdfplot(error_now);
        set(h, 'linewidth', 1);
        if area_id == 1
            title('S1');
        elseif area_id == 2
            title('S2');
        end % if area_id ==1
    end % for model_id = 1:2
    xticks([60, 100, 140, 180]);
    xlim([60, 200]);
    xlabel('minimum differences');
    ylabel('probability');
    hold off
    
end %area_id = 1:2
legend({'Sq(X)+Sq(Y)', 'Sq(X)+Sq(Y)+Sq(XY)'});
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width];
% Print
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'Sup_figure5_bp');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'Sup_figure5_bp', '-r300');
end



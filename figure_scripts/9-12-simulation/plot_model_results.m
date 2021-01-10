%% Plot model results
% model 
% aRect(X) + bSq(Y)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y)
% aRect(X) + bSq(Y) + cSq(XY)
% aRect(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
% aSq(X) + bSq(Y) + cSq(XY)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)

Rect_data = load('Results_Rect_top10.mat');
Sq_data = load('Results_Sq_top10.mat');

%% Rect histogram Best model id 
figure();clf
y_ = [];
yp_ = [];
for area_id = 1:2
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];    
    num_channels = size(data_each_channel, 2);
    % Histogram
    model_1 = sum(best_model_ids == 1);
    model_2 = sum(best_model_ids == 2);
    model_3 = sum(best_model_ids == 3);
    model_4 = sum(best_model_ids == 4);
    y_tmp = [model_1, model_2, model_3, model_4];
    y_ = [y_; y_tmp];
    % Probability
    yp_tmp = y_tmp ./ num_channels;
    yp_ = [yp_; yp_tmp];
end %area_id = 1:2
x_ = [1, 2, 3, 4];
bar(x_, yp_);    
xticks(x_);
xticklabels({'RX+RY', 'RX+RY+RXRY', 'RX+RY+RXY', 'all'});
%{
xticklabels({'Rect(X)+Rect(Y)', 'Rect(X)+Rect(Y)+Rect(X)Rect(Y)',...
        'Rect(X)+Rect(Y)+Rect(XY)',...
        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)+Rect(XY)'});
%}
ylim([0 1]);
ylabel('Probability');
sgtitle('Best Rect model');
legend({'S1', 'S2'});
set(gca, 'fontsize', 16);

%% Rect check each model
figure(); clf
for area_id = 1:2
    area = Rect_data.each_area(area_id).area;
    data_each_channel = Rect_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
    % + 
    hold on 
    subplot(1, 4, 1);
    h = histogram(errors(:,1));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 2);
    h = histogram(errors(:,2));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 3);
    h = histogram(errors(:,3));
    h.Normalization = 'probability';
    h.BinWidth = 1;    
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 4);
    h = histogram(errors(:,4));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');    
    ylabel('probability');
    hold off
end %area_id = 1:2

legend({'S1', 'S2'});


%% Sq histogram Best model id 
figure();clf
y_ = [];
yp_ = [];
for area_id = 1:2
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];    
    num_channels = size(data_each_channel, 2);
    model_1 = sum(best_model_ids == 1);
    model_2 = sum(best_model_ids == 2);
    model_3 = sum(best_model_ids == 3);
    model_4 = sum(best_model_ids == 4);
    y_tmp = [model_1, model_2, model_3, model_4];
    y_ = [y_; y_tmp];
    % Probability
    yp_tmp = y_tmp ./ num_channels;
    yp_ = [yp_; yp_tmp];
end %area_id = 1:2
x_ = [1, 2, 3, 4];
bar(x_, yp_);    
xticks(x_);
xticklabels({'SX+SY', 'SX+SY+SXRY', 'SX+SY+SXY', 'all'});
%{
xticklabels({'Sq(X)+Sq(Y)', 'Sq(X)+Sq(Y)+Sq(X)Sq(Y)',...
        'Sq(X)+Sq(Y)+Sq(XY)',...
        'Sq(X)+Sq(Y)+Sq(X)Sq(Y)+Sq(XY)'});
%}
ylim([0 1]);
ylabel('Probability');
sgtitle('Best Sq model');
legend({'S1', 'S2'});
set(gca, 'fontsize', 16);

%% Sqc check each model
figure(); clf
for area_id = 1:2
    area = Sq_data.each_area(area_id).area;
    data_each_channel = Sq_data.each_area(area_id).each_channel;
    best_model_ids = [data_each_channel.best_model_id];
    
    num_channels = size(data_each_channel, 2);
    errors = []; % channels x models
    for ch_id = 1:num_channels
        each_error = [data_each_channel(ch_id).each_model.error];
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
    % + 
    hold on 
    subplot(1, 4, 1);
    h = histogram(errors(:,1));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 2);
    h = histogram(errors(:,2));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 3);
    h = histogram(errors(:,3));
    h.Normalization = 'probability';
    h.BinWidth = 1;    
    xlabel('error');
    ylabel('probability');
    hold off
    
    hold on
    subplot(1, 4, 4);
    h = histogram(errors(:,4));
    h.Normalization = 'probability';
    h.BinWidth = 1;
    xlabel('error');    
    ylabel('probability');
    hold off
end %area_id = 1:2

legend({'S1', 'S2'});


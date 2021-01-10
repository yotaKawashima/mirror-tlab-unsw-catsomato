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
% aSq(X) + bSq(Y) + cSq(XY)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)

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
    
    % Store error data for each area
    error_each_area(area_id).errors = errors;
    
end %area_id = 1:2


% Two-sample Kolmogorov-Smirnov test for each model between S1 and S2
h = [];
p = [];
ks2stat = [];
for model_id = 1:4
    [h_, p_, ks2stat_] = kstest2(error_each_area(1).errors(:,model_id), ...
                                 error_each_area(2).errors(:,model_id));
    h = [h; h_];
    p = [p; p_];
    ks2stat = [ks2stat; ks2stat_];    
end % for model_id =1:4
ksstats_table = table(h, p, ks2stat);
disp(ksstats_table);

% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    errors = error_each_area(area_id).errors;
    
    for model_id = 1:4
        hold on
        subplot(2, 2, model_id);
        [h, stats] = cdfplot(errors(:,model_id));
        set(h, 'linewidth', 3);
        xlabel('minimum differences');
        ylabel('probability');
        hold off
        switch model_id
            case 1
                title_now_1 = 'Rect(X)+Rect(Y)'; 
                title_now_2 = '';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];
            case 2
                title_now_1 = 'Rect(X)+Rect(Y)';
                title_now_2 = '+Rect(X)Rect(Y)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];                
            case 3
                title_now_1 = 'Rect(X)+Rect(Y)'; 
                title_now_2 = '+Rect(XY)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];                
            case 4
                title_now_1 = 'Rect(X)+Rect(Y)'; 
                title_now_2 = '+Rect(X)Rect(Y)+Rect(XY)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];               
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];
        end % swithch model_id   
        title({title_now_1, title_now_2, title_now_3, title_now_4});
      
        xlim([40, 160]);
        set(gca, 'fontsize', 16);

    end % for model_id = 1:4

end % area_id = 1:2

legend({'S1', 'S2'});

%% Rect compare models for each area
% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    % Get error data
    errors = error_each_area(area_id).errors;
    
    subplot(2,1, area_id);
    hold on 
    for model_id = 1:4
        [h, stats] = cdfplot(errors(:,model_id));
        set(h, 'linewidth', 3);
        if area_id == 1
            title('S1');
            %set(h, 'color', [0, 0.4470, 0.7410]);
        elseif area_id == 2
            title('S2');
            %set(h, 'color', [0.8500, 0.3250, 0.0980]);
        end % if area_id ==1
        %{
        switch model_id
            case 1
                set(h, 'linestyle', '-');
                set(h, 'linewidth', 2.5);        
            case 2
                set(h, 'linestyle', ':');
                set(h, 'linewidth', 2.5);               
            case 3
                set(h, 'linestyle', ':');
                set(h, 'linewidth', 2.5);
                if area_id == 1
                    set(h, 'color', [0.4660, 0.6740, 0.1880]);
                elseif area_id == 2
                    set(h, 'color', [0.9290, 0.6940, 0.1250]);                    
                end % if area_id == 1                 
            case 4
                set(h, 'linestyle', '-');
                set(h, 'linewidth', 2.5);  
                if area_id == 1
                    set(h, 'color', [0.4660, 0.6740, 0.1880]);
                elseif area_id == 2
                    set(h, 'color', [0.9290, 0.6940, 0.1250]);                    
                end % if area_id == 1                 
        end % switch model_id  
        %}
    end % for model_id = 1:4
    xlim([40, 160]);
    xlabel('minimum differences');
    ylabel('probability');
    hold off
    legend({'Rect(X)+Rect(Y)', 'Rect(X)+Rect(Y)+Rect(X)Rect(Y)', ...
           'Rect(X)+Rect(Y)+Rect(XY)', ...
           'Rect(X)+Rect(Y)+Rect(X)Rect(Y)+Rect(XY)'});
    set(gca, 'fontsize', 16);
    
end %area_id = 1:2


%% Sq 
% aSq(X) + bSq(Y)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
% aSq(X) + bSq(Y) + cSq(XY)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
Sq_data = load('Results_Sq_top10.mat');

%% Sq compare S1 and S2 for each model
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
        errors = [errors; each_error];
    end % for ch_id = 1:num_channel
    
    % Store error data for each area
    error_each_area(area_id).errors = errors;
    
end %area_id = 1:2

% Two-sample Kolmogorov-Smirnov test for each model between S1 and S2
h = [];
p = [];
ks2stat = [];
for model_id = 1:4
    [h_, p_, ks2stat_] = kstest2(error_each_area(1).errors(:,model_id), ...
                                 error_each_area(2).errors(:,model_id));
    h = [h; h_];
    p = [p; p_];
    ks2stat = [ks2stat; ks2stat_];    
end % for model_id =1:4
ksstats_table = table(h, p, ks2stat);
disp(ksstats_table);

% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    errors = error_each_area(area_id).errors;
    
    for model_id = 1:4
        hold on
        subplot(2, 2, model_id);
        [h, stats] = cdfplot(errors(:,model_id));
        set(h, 'linewidth', 3);
        xlabel('minimum differences');
        ylabel('probability');
        hold off
        
        switch model_id
            case 1
                title_now_1 = 'Sq(X)+Sq(Y)'; 
                title_now_2 = '';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];   
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];                
            case 2
                title_now_1 = 'Sq(X)+Sq(Y)';
                title_now_2 = '+Sq(X)Sq(Y)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];                  
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];
            case 3
                title_now_1 = 'Sq(X)+Sq(Y)'; 
                title_now_2 = '+Sq(XY)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];                
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];
            case 4
                title_now_1 = 'Sq(X)+Sq(Y)'; 
                title_now_2 = '+Sq(X)Sq(Y)+Sq(XY)';
                title_now_3 = ['p = ', num2str(ksstats_table.p(model_id))];
                title_now_4 = ['ks = ', num2str(ksstats_table.ks2stat(model_id))];                
        end % swithch model_id
        title({title_now_1, title_now_2, title_now_3, title_now_4});
        
        xlim([60, 240]);
        set(gca, 'fontsize', 16);

    end % for model_id = 1:4

end % area_id = 1:2
legend({'S1', 'S2'});

%% Sq compare models for each area
% Plot culmtive probability distribution
figure(); clf
for area_id = 1:2
    % Get error data
    errors = error_each_area(area_id).errors;
    
    subplot(2,1, area_id);
    hold on 
    for model_id = 1:4
        [h, stats] = cdfplot(errors(:,model_id));
        set(h, 'linewidth', 3);
        if area_id == 1
            title('S1');
            %set(h, 'color', [0, 0.4470, 0.7410]);
        elseif area_id == 2
            title('S2');
            %set(h, 'color', [0.8500, 0.3250, 0.0980]);
        end % if area_id ==1
        %{
        switch model_id
            case 1
                set(h, 'linestyle', '-');
                set(h, 'linewidth', 2.5);        
            case 2
                set(h, 'linestyle', ':');
                set(h, 'linewidth', 2.5);               
            case 3
                set(h, 'linestyle', ':');
                set(h, 'linewidth', 2.5);
                if area_id == 1
                    set(h, 'color', [0.4660, 0.6740, 0.1880]);
                elseif area_id == 2
                    set(h, 'color', [0.9290, 0.6940, 0.1250]);                    
                end % if area_id == 1                 
            case 4
                set(h, 'linestyle', '-');
                set(h, 'linewidth', 2.5);
                if area_id == 1
                    set(h, 'color', [0.4660, 0.6740, 0.1880]);
                elseif area_id == 2
                    set(h, 'color', [0.9290, 0.6940, 0.1250]);                    
                end % if area_id == 1                 
                
        end % switch model_id  
        %}

    end % for model_id = 1:4
    xlim([60, 240]);
    xlabel('minimum differences');
    ylabel('probability');
    hold off
    legend({'Sq(X)+Sq(Y)', 'Sq(X)+Sq(Y)+Sq(X)Sq(Y)', ...
           'Sq(X)+Sq(Y)+Sq(XY)', ...
           'Sq(X)+Sq(Y)+Sq(X)Sq(Y)+Sq(XY)'});
    set(gca, 'fontsize', 16);
    
end %area_id = 1:2



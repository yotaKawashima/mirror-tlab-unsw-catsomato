%% Plot exemplar responses (Test)
cat_name = 'C20110808_R03';

% Load top 10 channels information.
load('Harmonics_IMs_top10.mat');

% Bipolar channel (Note this ID is not actual bipolar channel ID.)
bp_ch_id = 10;

img_fmt = '-dpdf';
fsize = 10; % font size
ftype = 'Arial'; % font type
x_width = 8; % fig width
y_width = 8; % fig height


foi_f1_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
f_harm_inter = sort([foi_f1_harm, foi_f2, foi_inter]);

%% Rectification 
% Load Rect data
Rect_data = load('Results_Rect_top10_v2.mat');

% Get error from S1
S1_channel_data = Rect_data.each_area(1).each_channel;
freqs = Rect_data.each_area(1).freqs;

% Find sessoin of our interest.
search_array = repmat(cat_name, size(S1_channel_data, 2), 1);
S1_session_array = data_top10_struct(1).session;
session_this_cat = prod(S1_session_array == search_array, 2);
session_ids_this_cat = find(session_this_cat == 1);
S1_channel_data_this_cat = S1_channel_data(1, session_ids_this_cat);
S1_bipolar_chs_this_cat = data_top10_struct(1).bipolar_ch(session_ids_this_cat);

disp('Rect');
disp(cat_name);
disp(['bipolar channel: ', num2str(S1_bipolar_chs_this_cat(bp_ch_id))]);

% Get logSNR data
recorded_logSNR = S1_channel_data_this_cat(bp_ch_id).logsnr;
model_1_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(1).model_logsnr;
model_2_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(2).model_logsnr;
model_3_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(3).model_logsnr;
model_4_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(4).model_logsnr;

% error 
error_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).error;
error_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).error;
error_3 = S1_channel_data_this_cat(bp_ch_id).each_model(3).error;
error_4 = S1_channel_data_this_cat(bp_ch_id).each_model(4).error;

% parameter 
params_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).parameters;
params_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).parameters;
params_3 = S1_channel_data_this_cat(bp_ch_id).each_model(3).parameters;
params_4 = S1_channel_data_this_cat(bp_ch_id).each_model(4).parameters;

% plot logSNR
figure();clf
subplot(5,1,1:3)
hold on
offset_val = 30;
% horizontal line
for i_hline = -2:2
    % plot horizontal line
    h = plot([0, 250], [i_hline, i_hline]*offset_val, 'k');
    % format vertical lines
    set(h, 'LineWidth', 0.5)
    % Remove legend for these lines
    h_annotation = get(h, 'Annotation');     
    set(get(h_annotation,'LegendInformation'),...
        'IconDisplayStyle','off');
    uistack(h,'bottom')
end % i_hline = -2:2
%Recorded data
plot(freqs, recorded_logSNR+offset_val*2, 'color', 'k', 'linewidth', 1);
%Rect(X)+Rect(Y)
plot(freqs, model_1_logSNR+offset_val, 'color', [89, 183, 250]/256, 'linewidth', 1); 
%Rect(X)+Rect(Y)+Rect(X)Rect(Y)
plot(freqs, model_2_logSNR-offset_val, 'color', [189, 89, 250]/256, 'linewidth', 1); 
%Rect(X)+Rect(Y)+Rect(XY)
plot(freqs, model_3_logSNR+0, 'color', [250, 125, 89]/256, 'linewidth', 1);     
%Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)
plot(freqs, model_4_logSNR-offset_val*2, 'color', [255, 195, 0]/256, 'linewidth', 1); 
% Add a line showing 20 
ind_line = plot([240, 240], [5, 25],'k', 'LineWidth', 1.5);
iline_annotation = get(ind_line, 'Annotation');     
set(get(iline_annotation,'LegendInformation'),...
    'IconDisplayStyle','off');
text(240, 15,'20');
xlim([0, 250]);
ylim([-20-2*offset_val, 20+2*offset_val]);
set(gca, 'YTick', offset_val*[-2 , -1, 0, 1, 2])
set(gca,'YTickLabel',[0, 0, 0, 0, 0]);
add_verticle_line(gca); 
xlabel('frequency [Hz]');
ylabel('logSNR [-]');


% Legend
%legend({'Observed';...
%        'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)'...
%        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';...
%        'Rect(X)+Rect(Y)+Rect(XY)';...
%        'Rect(X)+Rect(Y)';...
%        );
% title
title('Observed logSNR and model logSNR', 'fontsize', fsize);

hold off 

% Get points at f1 harmonics and intermodulation
[~, harm_inter_inds] = find_closest(freqs, f_harm_inter);  
x_points = freqs(harm_inter_inds);
diff_1 = recorded_logSNR - model_1_logSNR;
diff_3 = recorded_logSNR - model_3_logSNR;
diff_2 = recorded_logSNR - model_2_logSNR;
diff_4 = recorded_logSNR - model_4_logSNR;

subplot(5,1,5)    
hold on 
%Difference Recorded data
%Rect(X)+Rect(Y)
scatter(x_points-1.25, diff_1(harm_inter_inds), 30, [89, 183, 250]/256, 'filled'); 
%Rect(X)+Rect(Y)+Rect(XY)
scatter(x_points-0.625, diff_3(harm_inter_inds), 30, [250, 125, 89]/256, 'filled');     
%Rect(X)+Rect(Y)+Rect(X)Rect(Y)
scatter(x_points+0.625, diff_2(harm_inter_inds), 30, [189, 89, 250]/256, 'filled'); 
%Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)
scatter(x_points+1.25, diff_4(harm_inter_inds), 30, [255, 195, 0]/256, 'filled'); 
xlim([0,250]);
ylim([-25, 25]);
add_verticle_line(gca);
xlabel('frequency [Hz]');
ylabel('difference [-]');

% Legend
legend({'Rect(X)+Rect(Y)';...
        'Rect(X)+Rect(Y)+Rect(XY)';...
        'Rect(X)+Rect(Y)+Rect(X)Rect(Y)';...
        'Rect(X)+Rect(Y)+Rect(XY)+Rect(X)Rect(Y)'});
% title
title('Difference:(observed logSNR) - (model logSNR)');
hold off
set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');    
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [5, 5, x_width, y_width*2.5];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'figure10');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'figure10', '-r300');
end

%% Square
% Load Sq data
Sq_data = load('Results_Sq_top10_v2.mat');

% Get error from S1
S1_channel_data = Sq_data.each_area(1).each_channel;
freqs = Sq_data.each_area(1).freqs;

% Find sessoin of our interest.
search_array = repmat(cat_name, size(S1_channel_data, 2), 1);
S1_session_array = data_top10_struct(1).session;
session_this_cat = prod(S1_session_array == search_array, 2);
session_ids_this_cat = find(session_this_cat == 1);
S1_channel_data_this_cat = S1_channel_data(1, session_ids_this_cat);
S1_bipolar_chs_this_cat = data_top10_struct(1).bipolar_ch(session_ids_this_cat);

disp('Sq');
disp(cat_name);
disp(['bipolar channel: ', num2str(S1_bipolar_chs_this_cat(bp_ch_id))]);

% Get logSNR data
recorded_logSNR = S1_channel_data_this_cat(bp_ch_id).logsnr;
model_1_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(1).model_logsnr;
model_2_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(2).model_logsnr;

% error 
error_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).error;
error_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).error;

% parameter 
params_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).parameters;
params_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).parameters;

% plot logSNR
figure();clf

subplot(4,1,1:2)  
hold on
offset_val = 30;
% horizontal line
for i_hline = -1:1
    % plot horizontal line
    h = plot([0, 250], [i_hline, i_hline]*offset_val, 'k');
    % format vertical lines
    set(h, 'LineWidth', 0.5)
    % Remove legend for these lines
    h_annotation = get(h, 'Annotation');     
    set(get(h_annotation,'LegendInformation'),...
        'IconDisplayStyle','off');
    uistack(h,'bottom')
end % i_hline = -1:1
%Recorded data
plot(freqs, recorded_logSNR+offset_val, 'color', 'k', 'linewidth', 1);
%Sq(X)+Sq(Y)
plot(freqs, model_1_logSNR+0, 'color', [89, 183, 250]/256, 'linewidth', 1); 
%Sq(X)+Sq(Y)+Sq(XY)
plot(freqs, model_2_logSNR-offset_val, 'color', [250, 125, 89]/256, 'linewidth', 1);     
% Add a line showing 20 
ind_line = plot([230, 230], [5, 25],'k', 'LineWidth', 1.5);
iline_annotation = get(ind_line, 'Annotation');     
set(get(iline_annotation,'LegendInformation'),...
    'IconDisplayStyle','off');
text(235, 15,'20');
xlim([0, 250]);
ylim([-20-offset_val, 20+offset_val]);
set(gca, 'YTick', offset_val*[ -1, 0, 1])
set(gca,'YTickLabel',[0, 0, 0]);
add_verticle_line(gca); 
xlabel('frequency [Hz]');
ylabel('logSNR [-]');

% Legend
%legend({'Observed'; 'Sq(X)+Sq(Y)+Sq(XY)'; 'Sq(X)+Sq(Y)'});

% title
title('Observed logSNR and model logSNR');
hold off 

% Get points at f1 harmonics and intermodulation
[~, harm_inter_inds] = find_closest(freqs, f_harm_inter);  
x_points = freqs(harm_inter_inds);
diff_1 = recorded_logSNR - model_1_logSNR;
diff_2 = recorded_logSNR - model_2_logSNR;

subplot(4,1,4)  
hold on
%Sq(X)+Sq(Y)
scatter(x_points-0.75, diff_1(harm_inter_inds), 30, [89, 183, 250]/256, 'filled'); 
%Sq(X)+Sq(Y)+Sq(XY)
scatter(x_points+0.75, diff_2(harm_inter_inds), 30, [250, 125, 89]/256, 'filled');     
xlim([0,250]);
ylim([-25, 25]);
add_verticle_line(gca); 
xlabel('frequency [Hz]');
ylabel('difference [-]');

% Legend
legend({'Sq(X)+Sq(Y)'; 'Sq(X)+Sq(Y)+Sq(XY)'});

% title
title('Difference:(observed logSNR) - (model logSNR)');
hold off 

set(findall(gcf,'-property','FontSize'), 'FontSize', fsize);
set(findall(gcf,'-property','FontName'), 'FontName', ftype);
%set(gcf,'color','w');    
set(gcf,'renderer','Painters');
f=gcf;
f.Units = 'centimeters';
f.Position = [10, 10, x_width, y_width*2];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(gcf, img_fmt, 'Sup_figure4');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'Sup_figure4', '-r300');
end

%% functions
function add_verticle_line(gca)
%Add verticle lines indicating fundamentals, harmonics, and IMs.
%  input gca
    % Frequencies of interest
    foi_f1_and_harm = 23:23:250;
    foi_f2 = 200;
    foi_inter = sort([177:-23:0 200+23:23:250]);
    foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

    
    set(gca, 'fontsize', 16);
    
    % plot vertical lines at frequencies of interest
    % grab axis variables
    ylims = get(gca, 'YLim');
    v_ylims = [ylims(1), (ylims(2)-ylims(1))*0.05 + ylims(1)];
    v1 = numel(foi_f1_and_harm);
    v2 = numel(foi_inter);
    v = zeros(1, numel(foi));
    % plot vertical lines
    v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
    v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
    v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
    % format vertical lines
    set(v, 'LineWidth', 1.5 )
    % Remove legend for these lines
    v_annotation = get(v, 'Annotation');
    for i_vline = 1:length(v)
        set(get(v_annotation{i_vline,1},'LegendInformation'),...
            'IconDisplayStyle','off');
    end % i_vline = 1:length(v)
    
    hold off;    
end

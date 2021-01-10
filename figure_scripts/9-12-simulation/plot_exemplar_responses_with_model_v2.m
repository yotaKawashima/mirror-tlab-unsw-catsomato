%% Plot exemplar responses (Test)
cat_name = 'C20110808_R03';

% Load top 10 channels information.
load('Harmonics_IMs_top10.mat');

% Bipolar channel (Note this ID is not actual bipolar channel ID.)
bp_ch_id = 10;

%% Rectification 
% Load Rect data
Rect_data = load('Results_Rect_top10.mat');

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

%Rect(X)+Rect(Y)
subplot(4,1,1)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_1_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2); 
add_verticle_line(gca);    
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;
% title
title(['error:', num2str(error_1), ', parameters: (a,b)=', ...
       '(', num2str(params_1(1)), ', ', num2str(params_1(2)), ')']);
hold off 

%Rect(X)+Rect(Y)+Rect(XY)
subplot(4,1,2)    
hold on 
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_3_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2);     
add_verticle_line(gca);
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;
% title
title(['error:', num2str(error_3), ', parameters: (a,b,c)=', ...
       '(', num2str(params_3(1)), ', ', num2str(params_3(2)), ...
       ', ', num2str(params_3(3)), ')']);
hold off
    
%Rect(X)+Rect(Y)+Rect(X)Rect(Y)
subplot(4,1,3)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);    
plot(freqs, model_2_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2); 
add_verticle_line(gca);
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;    
% title
title(['error:', num2str(error_2), ', parameters: (a,b,c)=', ...
       '(', num2str(params_2(1)), ', ', num2str(params_2(2)), ...
       ', ', num2str(params_2(3)), ')']);
hold off

%Rect(X)+Rect(Y)+Rect(X)Rect(Y)+Rect(XY)
subplot(4,1,4)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_4_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2);
add_verticle_line(gca);  
% Legend
legend({'Observed data'; 'Model'});
axis tight;
% title
title(['error:', num2str(error_4), ', parameters: (a,b,c,d)=', ...
       '(', num2str(params_4(1)), ', ', num2str(params_4(2)), ...
       ', ', num2str(params_4(3)), ', ', num2str(params_4(4)), ')']);   
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
hold off

%% Square
% Load Sq data
Sq_data = load('Results_Sq_top10.mat');

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

%Sq(X)+Sq(Y)
subplot(4,1,1)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_1_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2); 
add_verticle_line(gca);    
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;
% title
title(['error:', num2str(error_1), ', parameters: (a,b)=', ...
       '(', num2str(params_1(1)), ', ', num2str(params_1(2)), ')']);
hold off 

%Sq(X)+Sq(Y)+Sq(XY)
subplot(4,1,2)    
hold on 
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_3_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2);     
add_verticle_line(gca);
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;
% title
title(['error:', num2str(error_3), ', parameters: (a,b,c)=', ...
       '(', num2str(params_3(1)), ', ', num2str(params_3(2)), ...
       ', ', num2str(params_3(3)), ')']);
hold off
    
%Sq(X)+Sq(Y)+Sq(X)Sq(Y)
subplot(4,1,3)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);    
plot(freqs, model_2_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2); 
add_verticle_line(gca);
% Legend
legend({'Observed data'; 'Model'});
set(gca,'xticklabel',[]);
axis tight;    
% title
title(['error:', num2str(error_2), ', parameters: (a,b,c)=', ...
       '(', num2str(params_2(1)), ', ', num2str(params_2(2)), ...
       ', ', num2str(params_2(3)), ')']);
hold off

%Sq(X)+Sq(Y)+Sq(X)Sq(Y)+Rect(XY)
subplot(4,1,4)
hold on
%Recorded data
plot(freqs, recorded_logSNR, 'color', 'k', 'linewidth', 3);
plot(freqs, model_4_logSNR, 'color', [204, 132, 63]/256, 'linewidth', 2);
add_verticle_line(gca);  
% Legend
legend({'Observed data'; 'Model'});
axis tight;
% title
title(['error:', num2str(error_4), ', parameters: (a,b,c,d)=', ...
       '(', num2str(params_4(1)), ', ', num2str(params_4(2)), ...
       ', ', num2str(params_4(3)), ', ', num2str(params_4(4)), ')']);
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
hold off


%% functions
function add_verticle_line(gca)
%Add verticle lines indicating fundamentals, harmonics, and IMs.
%  input gca
    % Frequencies of interest (Just for visualisatoin)
    foi_f1_and_harm = 23:23:250;
    foi_f2 = 200;
    foi_inter = sort([177:-23:0 200+23:23:250]);
    foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

    
    set(gca, 'fontsize', 16);
    
    % plot vertical lines at frequencies of interest
    % grab axis variables
    ylims = get(gca, 'YLim');
    v_ylims = [-15, -10];
    v1 = numel(foi_f1_and_harm);
    v2 = numel(foi_inter);
    v = zeros(1, numel(foi));
    % plot vertical lines
    v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
    v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
    v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
    % format vertical lines
    set(v, 'LineWidth', 2, 'LineStyle', '--')
    % Remove legend for these lines
    v_annotation = get(v, 'Annotation');
    for i_vline = 1:length(v)
        set(get(v_annotation{i_vline,1},'LegendInformation'),...
            'IconDisplayStyle','off');
    end % i_vline = 1:length(v)
    
    set(gca, 'XLim', [0,250]);
    set(gca, 'YLim', [-15, 35]);
    
    hold off;    
end

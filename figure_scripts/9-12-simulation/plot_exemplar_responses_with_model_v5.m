%% Set overall variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

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
%Rect(X)+Rect(Y)+Rect(XY)
plot(freqs, model_3_logSNR+0, 'color', [250, 125, 89]/256, 'linewidth', 1);     
%Rect(X)+Rect(Y)+Rect(X)Rect(Y)
plot(freqs, model_2_logSNR-offset_val, 'color', [189, 89, 250]/256, 'linewidth', 1); 
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
ylabel('logSNR [dB]');


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
ylabel('difference [dB]');

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

%% Half Squaring
% Load HSq data
HSq_data = load('Results_HSq_top10_v2.mat');

% Get error from S1
S1_channel_data = HSq_data.each_area(1).each_channel;
freqs = HSq_data.each_area(1).freqs;

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
%HSq(X)+HSq(Y)
plot(freqs, model_1_logSNR+offset_val, 'color', [89, 183, 250]/256, 'linewidth', 1); 
%HSq(X)+HSq(Y)+HSq(XY)
plot(freqs, model_3_logSNR+0, 'color', [250, 125, 89]/256, 'linewidth', 1);     
%HSq(X)+HSq(Y)+HSq(X)HSq(Y)
plot(freqs, model_2_logSNR-offset_val, 'color', [189, 89, 250]/256, 'linewidth', 1); 
%HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)
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
ylabel('logSNR [dB]');


% Legend
%legend({'Observed';...
%        'HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)'...
%        'HSq(X)+HSq(Y)+HSq(X)HSq(Y)';...
%        'HSq(X)+HSq(Y)+HSq(XY)';...
%        'HSq(X)+HSq(Y)';...
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
%HSq(X)+HSq(Y)
scatter(x_points-1.25, diff_1(harm_inter_inds), 30, [89, 183, 250]/256, 'filled'); 
%HSq(X)+HSq(Y)+HSq(XY)
scatter(x_points-0.625, diff_3(harm_inter_inds), 30, [250, 125, 89]/256, 'filled');     
%HSq(X)+HSq(Y)+HSq(X)HSq(Y)
scatter(x_points+0.625, diff_2(harm_inter_inds), 30, [189, 89, 250]/256, 'filled'); 
%HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)
scatter(x_points+1.25, diff_4(harm_inter_inds), 30, [255, 195, 0]/256, 'filled'); 
xlim([0,250]);
ylim([-25, 25]);
add_verticle_line(gca);
xlabel('frequency [Hz]');
ylabel('difference [dB]');

% Legend
legend({'HSq(X)+HSq(Y)';...
        'HSq(X)+HSq(Y)+HSq(XY)';...
        'HSq(X)+HSq(Y)+HSq(X)HSq(Y)';...
        'HSq(X)+HSq(Y)+HSq(XY)+HSq(X)HSq(Y)'});
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
    print(gcf, img_fmt, 'Sup_figure4');
elseif img_fmt == "-dtiff"
    print(gcf, img_fmt, 'Sup_figure4', '-r300');
end

%% Check model parameters
%{
% Sampling rate [Hz]
sampling_rate = 10000; 
model_4_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(4).model_logsnr;
model_4_params = S1_channel_data_this_cat(bp_ch_id).each_model(4).parameters;
%model_4_signal_from_params = ...
%    model_all(model_4_params(1), model_4_params(2),...
%              model_4_params(3), model_4_params(4));
model_4_signal_from_params = ...
    model_all(-0.85, 35,...
              1.3, 1.4e+308);

[model_4_logSNR_from_params, ~, freqs]...
    = compute_logsnrs_y(model_4_signal_from_params, sampling_rate);
                              
figure();clf
hold on;
plot(freqs, model_4_logSNR);
plot(freqs, model_4_logSNR_from_params);
hold off;
%}

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
%{
function model = model_plus(a, b)
% aHSq(X) + bHsq(Y)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Half Square (x^2 if x >0 otherwise 0)
    HSq_comp1 = zeros(size(comp1));
    HSq_comp1(comp1 > 0) = comp1(comp1 > 0).^2;
    HSq_comp1(comp1 <=0 ) = 0;
    
    HSq_comp2 = zeros(size(comp2));
    HSq_comp2(comp2 > 0) = comp2(comp2 > 0).^2;
    HSq_comp2(comp2 <=0 ) = 0;    

    model = (a*HSq_comp1) + (b*HSq_comp2);
end

function model = model_without_HSqXY(a, b, c)
% aHSq(X) + bHSq(Y) + cHSq(X)HSq(Y)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Half Square (x^2 if x >0 otherwise 0)
    HSq_comp1 = zeros(size(comp1));
    HSq_comp1(comp1 > 0) = comp1(comp1 > 0).^2;
    HSq_comp1(comp1 <=0 ) = 0;
    
    HSq_comp2 = zeros(size(comp2));
    HSq_comp2(comp2 > 0) = comp2(comp2 > 0).^2;
    HSq_comp2(comp2 <=0 ) = 0;    
        
    % model 
    model = a*HSq_comp1 + b*HSq_comp2 + c*HSq_comp1 .* HSq_comp2;
end

function model = model_without_HSqXHSqY(a, b, c)
% aHSq(X) + bHSq(Y) + cHSq(XY)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Half Square (x^2 if x >0 otherwise 0)
    HSq_comp1 = zeros(size(comp1));
    HSq_comp1(comp1 > 0) = comp1(comp1 > 0).^2;
    HSq_comp1(comp1 <=0 ) = 0;
    
    HSq_comp2 = zeros(size(comp2));
    HSq_comp2(comp2 > 0) = comp2(comp2 > 0).^2;
    HSq_comp2(comp2 <=0 ) = 0;  
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    HSq_mult_sines = zeros(size(mult_sines));
    HSq_mult_sines(mult_sines > 0) = ...
        HSq_mult_sines(mult_sines > 0).^2;
    HSq_mult_sines(mult_sines <= 0) = 0;
    
    % model 
    model = a*HSq_comp1 + b*HSq_comp2 + c*HSq_mult_sines;

end

function model = model_all(a, b, c, d)
% aHSq(X) + bHSq(Y) + cHSq(X)HSq(Y) + dHSq(XY)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);
    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);
    
    % Half Square (x^2 if x >0 otherwise 0)
    HSq_comp1 = zeros(size(comp1));
    HSq_comp1(comp1 > 0) = comp1(comp1 > 0).^2;
    HSq_comp1(comp1 <=0 ) = 0;
    
    HSq_comp2 = zeros(size(comp2));
    HSq_comp2(comp2 > 0) = comp2(comp2 > 0).^2;
    HSq_comp2(comp2 <=0 ) = 0;  
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    HSq_mult_sines = zeros(size(mult_sines));
    HSq_mult_sines(mult_sines > 0) = ...
        HSq_mult_sines(mult_sines > 0).^2;
    HSq_mult_sines(mult_sines <= 0) = 0;
    
    % model 
    model = a*HSq_comp1 + b*HSq_comp2 + c*HSq_comp1 .* HSq_comp2 ...
            + d*HSq_mult_sines;
end

function [logsnrs, powers, freqs] = compute_logsnrs_y(data, sampling_rate)
% Compute power spectrum and return logSNR.
% The same logsnr computation as Renee's snrtosurrounds.m, but this code is
% more easy to follow what we are doing.
%
%   The neighbourhood region is illustrated in the following diagram.
%   F=frequency of interest
%   A1=F - outerdiff       F - 3 Hz
%   B1=F - innerdiff       F - 1 Hz
%   C1=F - halfbandwidth   F - 0.5Hz
%   C2=F + halfbandwidth   F + 0.5Hz
%   B2=F + innerdiff       F + 1 Hz
%   A2=F + outerdiff       F + 3 Hz
%   
%   L = -1/(# of samples in low region)
%   H = +1/(# of samples in high region)
%   
%      LLLLLLLLLLL000HHHHH000LLLLLLLLLLL    <- Mask 
%   ... ......... ... ... ... ... ......... ..  <- data points 
%      |         |   |   |   |   |         |
%      A1        B1 C1   F   C2  B2        A2   <-freq
%      |         |   |   |   |   |         |
%     -3        -1  -0.5 0 +0.5  +1        +3
%                     ___ ___              
%    __           ___|       |___           __ 
%      |_________|               |_________|
%          low         high          low
%       
%   (This is just visualisation. # of points is not the same
%   as real data.

% Written by Yota. 
    % Parameters given to chronux mtspectrumc.m.
    params.Fs = sampling_rate;
    params.tapers = [1 1]; % hbw = (k+1)/2T  k=# of tapers, T=length of time
    params.pad = 1;   
    % Here, k = 1, T = 2s. Then, hbw = (1+1)/2*2 = 0.5Hz 
    halfbandwidth = 0.5; 
    
    % Compute power spectrum
    [p_part, f_part] = mtspectrumc(data, params);
    powers = 10*log10(p_part');
    freqs = f_part';

    % Parameter for computing logSNR
    innerdiff = 1;
    outerdiff = 3;
    half_width_low = outerdiff - innerdiff;
    maxf = 250;

    % create the mask
    % first, find the number of indicies I need to be 1 Hz
    % assume they're all evenly spaced. 
    spacing = freqs(2) - freqs(1);
    nSamp = round(1/spacing);    % # of indicies within 1Hz width. 
    nSamp_hbw = floor(nSamp*halfbandwidth); % # of indicies in hbw
    half_nSamp_low = half_width_low*nSamp; % # of indicies within low range
    
    % mask covering from A1 to A2.
    mLen = 2*(2*nSamp_hbw+half_nSamp_low)+1; % This could be wrong
    mask = zeros(1, mLen);

    % neighbour region outside C1-F-C2
    % -1/(# of Samples in neighbor region) 
    nSamp_low = 2*half_nSamp_low;      % # of indices in low
    lomask = 1/nSamp_low;                    % -1/(# of samples in low)
    mask(1:half_nSamp_low) = -lomask;         % from A1 to B1
    mask(end-half_nSamp_low+1:end) = -lomask; % from B2 to A2

    % central region C1-F-C2 (i.e. tagged freq)
    % # of samples in high 
    nSamp_high = nSamp_hbw*2 + 1; % Here, 1 is for data point of F
    himask = 1/nSamp_high;        %  +1/(# of samples in high) 
    m = ceil(mLen/2);             % Get the center of data point.
    mask(m-nSamp_hbw:m+nSamp_hbw) = himask;
    
    % reduce to only frequencies of interest
    [~, maxf_ind] = find_closest(freqs, maxf);
    freqs_interest = freqs(1:maxf_ind);
    
    % first frequency index at hbw*nS 
    % Some data points in the beginning and the end should be removed
    % because neighboring does not exist for them.
    % We can do zero-padding (to signal not to mask) to solve this. 
    % first frequency index at 3*nSample. 
    freqs = freqs_interest(nSamp_hbw*nSamp:(end-nSamp_hbw*nSamp+1));
    nReps = length(freqs);
    
    logsnrs = zeros(1, nReps);

    for index = 1:nReps
        % zero pad mask
        maskm = [zeros(1, index-1) mask zeros(1,maxf_ind-mLen-index+1)];
        logsnrs(1, index) = dot(maskm, powers(1:maxf_ind));
    end
    
    % Extract power within a range for 0Hz to 250Hz
    powers_interest = powers(1:maxf_ind);
    powers = powers_interest(nSamp_hbw*nSamp:(end-nSamp_hbw*nSamp+1));
    
end
%}

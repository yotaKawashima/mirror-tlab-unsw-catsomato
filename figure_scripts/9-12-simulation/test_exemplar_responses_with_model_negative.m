%% Plot exemplar responses (Test)
% check whether chaning sing of coefficient (from positive to negative)
% results in the same or different outputs.

sampling_rate = 10000; 

cat_name = 'C20110808_R03';

% Load top 10 channels information.
load('Harmonics_IMs_top10.mat');

% Bipolar channel (Note this ID is not actual bipolar channel ID.)
bp_ch_id = 10;

% Frequencies of interest (harmonics and intermodulaitons)
foi_f1_harm = 46:23:250;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi_now = sort([foi_f1_harm, foi_inter]);

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

% Get logSNR data with positive coefficients
positive_model_1_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(1).model_logsnr;
positive_model_2_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(2).model_logsnr;
positive_model_3_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(3).model_logsnr;
positive_model_4_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(4).model_logsnr;

% parameter 
positive_params_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).parameters;
positive_params_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).parameters;
positive_params_3 = S1_channel_data_this_cat(bp_ch_id).each_model(3).parameters;
positive_params_4 = S1_channel_data_this_cat(bp_ch_id).each_model(4).parameters;

% Flip parameter = 
negative_params_1 = -positive_params_1;
negative_params_2 = -positive_params_2;
negative_params_3 = -positive_params_3;
negative_params_4 = -positive_params_4;

% logSNR with negative coefficeints 
negative_signal_1 = model_plus(negative_params_1(1), negative_params_1(2));
negative_signal_2 = model_without_RXY(negative_params_2(1),...
                                      negative_params_2(2),... 
                                      negative_params_2(3));
negative_signal_3 = model_without_RXRY(negative_params_3(1),...
                                       negative_params_3(2),... 
                                       negative_params_3(3));
negative_signal_4 = model_all(negative_params_4(1), negative_params_4(2),... 
                              negative_params_4(3), negative_params_4(4));

 % Compute logSNR for the best model
[negative_model_1_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_1, ...
                                                    sampling_rate);
[negative_model_2_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_2, ...
                                                    sampling_rate);
[negative_model_3_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_3, ...
                                                    sampling_rate);
[negative_model_4_logSNR, ~, freqs_] = compute_logsnrs_y(negative_signal_4, ...
                                                    sampling_rate);                                                
        
% Compare negative vs positive.
diff_model_1 = negative_model_1_logSNR - positive_model_1_logSNR;
diff_model_2 = negative_model_2_logSNR - positive_model_2_logSNR;
diff_model_3 = negative_model_3_logSNR - positive_model_3_logSNR;
diff_model_4 = negative_model_4_logSNR - positive_model_4_logSNR;

figure();
hold on;
plot(freqs, diff_model_1);
plot(freqs, diff_model_2);
plot(freqs, diff_model_3);
plot(freqs, diff_model_4);
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
positive_model_1_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(1).model_logsnr;
positive_model_2_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(2).model_logsnr;
positive_model_3_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(3).model_logsnr;
positive_model_4_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(4).model_logsnr;

% parameter 
positive_params_1 = S1_channel_data_this_cat(bp_ch_id).each_model(1).parameters;
positive_params_2 = S1_channel_data_this_cat(bp_ch_id).each_model(2).parameters;
positive_params_3 = S1_channel_data_this_cat(bp_ch_id).each_model(3).parameters;
positive_params_4 = S1_channel_data_this_cat(bp_ch_id).each_model(4).parameters;

% Flip parameter = 
negative_params_1 = -positive_params_1;
negative_params_2 = -positive_params_2;
negative_params_3 = -positive_params_3;
negative_params_4 = -positive_params_4;

% logSNR with negative coefficeints 
negative_signal_1 = model_plus_Sq(negative_params_1(1), negative_params_1(2));
negative_signal_2 = model_without_RXY_Sq(negative_params_2(1),...
                                      negative_params_2(2),... 
                                      negative_params_2(3));
negative_signal_3 = model_without_RXRY_Sq(negative_params_3(1),...
                                       negative_params_3(2),... 
                                       negative_params_3(3));
negative_signal_4 = model_all_Sq(negative_params_4(1), negative_params_4(2),... 
                              negative_params_4(3), negative_params_4(4));

 % Compute logSNR for the best model
[negative_model_1_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_1, ...
                                                    sampling_rate);
[negative_model_2_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_2, ...
                                                    sampling_rate);
[negative_model_3_logSNR, ~, ~] = compute_logsnrs_y(negative_signal_3, ...
                                                    sampling_rate);
[negative_model_4_logSNR, ~, freqs_] = compute_logsnrs_y(negative_signal_4, ...
                                                    sampling_rate);                                                
        
% Compare negative vs positive.
diff_model_1 = negative_model_1_logSNR - positive_model_1_logSNR;
diff_model_2 = negative_model_2_logSNR - positive_model_2_logSNR;
diff_model_3 = negative_model_3_logSNR - positive_model_3_logSNR;
diff_model_4 = negative_model_4_logSNR - positive_model_4_logSNR;


figure();
hold on;
plot(freqs, diff_model_1);
plot(freqs, diff_model_2);
plot(freqs, diff_model_3);
plot(freqs, diff_model_4);
hold off

%% Functions
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
    
    hold off;    
end


%% Functions for computing logSNR
function error = compute_error(model_signal, datalogsnrs, sampling_rate)
% Compute error
%   Input: sampling_rate
%   Output: error

    % Frequencies of interest (Just for visualisatoin)
    foi_f1_and_harm = 23:23:250;
    foi_f2 = 200;
    foi_inter = sort([177:-23:0 200+23:23:250]);
    foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);            

    % Compute logSNR
    [modellogsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_signal, sampling_rate);
    
    [~, foi_inds] = find_closest(freqs, foi);
    error = sum(abs(modellogsnrs(1, foi_inds) - datalogsnrs(1, foi_inds)));
    
    % Check how log power looks like
    %{
    figure();clf
    [~, maxf_ind] = find_closest(freqs, 250);
    plot(freqs(1:maxf_ind), powers_tmp(1:maxf_ind))
    close gcf;    
    %}
    
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

%% Model 
function model = model_plus(a, b)
% aRect(X) + bRect(Y)
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

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = (a*Rect_comp1) + (b*Rect_comp2);
end

function model = model_without_RXY(a, b, c)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y)
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

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
        
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_comp1 .* Rect_comp2;
end

function model = model_without_RXRY(a, b, c)
% aRect(X) + bRect(Y) + cRect(XY)
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

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Rect_mult_sines = mult_sines;
    Rect_mult_sines(mult_sines <=0 ) = 0;
    
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_mult_sines;

end

function model = model_all(a, b, c, d)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY)
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

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Rect_mult_sines = mult_sines;
    Rect_mult_sines(mult_sines <=0 ) = 0;
    
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_comp1 .* Rect_comp2 ...
            + d*Rect_mult_sines;
end


%% Model 
function model = model_plus_Sq(a, b)
% aSq(X) + bSq(Y)
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
    
    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;
    
    % model
    model = (a*Sq_comp1) + (b*Sq_comp2);
end

function model = model_without_RXY_Sq(a, b, c)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
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

    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;
        
    % model 
    model = a*Sq_comp1 + b*Sq_comp2 + c*Sq_comp1 .* Sq_comp2;
end

function model = model_without_RXRY_Sq(a, b, c)
% aSq(X) + bSq(Y) + cSq(XY)
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

    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;

    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Sq_mult_sines(sin(mult_sines) > 0) = 1;
    Sq_mult_sines(sin(mult_sines) <=0) = -1;
    
    % model 
    model = a*Sq_comp1 + b*Sq_comp2 + c*Sq_mult_sines;

end

function model = model_all_Sq(a, b, c, d)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y) + dSq(XY)
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

    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;

    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Sq_mult_sines(sin(mult_sines) > 0) = 1;
    Sq_mult_sines(sin(mult_sines) <=0) = -1;
    
    % model 
    model = a*Sq_comp1 + b*Sq_comp2 + c*Sq_comp1 .* Sq_comp2 ...
            + d*Sq_mult_sines;
end

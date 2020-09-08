%% Plot exemplar responses (Test)
cat_name = 'C20110808_R03';


% Load top 10 channels information.
load('Harmonics_IMs_top10.mat');

%% Rectification 
% Load Rect data
Rect_data = load('Results_Rect_top10.mat');

% Get error from each area
S1_channel_data = Rect_data.each_area(1).each_channel;
freqs = Rect_data.each_area(1).freqs;

% Find sessoin of our interest.
search_array = repmat(cat_name, size(S1_channel_data, 2), 1);
S1_session_array = data_top10_struct(1).session;
session_this_cat = prod(S1_session_array == search_array, 2);
session_ids_this_cat = find(session_this_cat == 1);
S1_channel_data_this_cat = S1_channel_data(1, session_ids_this_cat);
S1_bipolar_chs_this_cat = data_top10_struct(1).bipolar_ch(session_ids_this_cat);

%%
for id =10:10
    disp(cat_name);
    disp(['bipolar channel: ', num2str(S1_bipolar_chs_this_cat(id))]);

    % Get logSNR data
    recorded_logSNR = S1_channel_data_this_cat(id).logsnr;
    model_1_logSNR = S1_channel_data_this_cat(id).each_model(1).model_logsnr;
    model_2_logSNR = S1_channel_data_this_cat(id).each_model(2).model_logsnr;
    model_3_logSNR = S1_channel_data_this_cat(id).each_model(3).model_logsnr;
    model_4_logSNR = S1_channel_data_this_cat(id).each_model(4).model_logsnr;


    % plot logSNR
    figure();clf
    hold on
    plot(freqs, recorded_logSNR, 'k'); %Recorded data
    %plot(freqs, model_1_logSNR); %Rect(X)+Rect(Y)
    %plot(freqs, model_2_logSNR); %Rect(X)+Rect(Y)+Rect(X)Rect(Y)
    %plot(freqs, model_3_logSNR); %Rect(X)+Rect(Y)+Rect(XY)
    plot(freqs, model_4_logSNR); %Rect(X)+Rect(Y)+Rect(X)Rect(Y)+Rect(XY)
    hold off;
end
%% Square

signal = model_all(0.0136316951162585,0.5,0.0318551123612303,3.32187397540912e-07);

% Sampling rate [Hz]
sampling_rate = 10000; 
[logsnrs, powers, freqs] = compute_logsnrs_y(signal, sampling_rate); 
hold on
plot(freqs, logsnrs);
hold off


%% Functions
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

%% Models
function model = model_all(a, b, c, d)
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




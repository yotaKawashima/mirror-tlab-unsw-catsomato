%% Fitting some Model to actual data for comparison. (Test version 2)
% Here, we focus on only multiplication model, and see how the singal
% itself look like with different parameters.

%% Set variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));


catname = 'C20110808_R03';
area = 'S1';
% If you included unipolar channels in your stored data. You need to be
% careful that bipolar channel id starts at 1011 and at 65 for S1 and S2.
bp_chan_1 = 256; % bipolar channel
bp_chan_2 = 152; % bipolar channel

% Frequencies of interest (Just for visualisatoin)
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

% frequency
f1 = 23;
f2 = 200;

% Sampling rate [Hz]
sampling_rate = 10000; 

% time 0s to 2s.
time_s = 0; % Start [s]
time_e = 2; % End [s]
times = time_s:1/sampling_rate:time_e;

% parameter for model
a_range = [0:0.01:0.1];
b_range = [0:0.1:5];

%% Model
%{
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
%}
%model = @(a,b)a*Rect_comp1 + b*Rect_comp2;
%model = @(a,b)(a*Rect_comp1) .* (b*Rect_comp2);
%model = @(a, b)model_plus(a,b); % a and b are inslide the Rect()
model = @(a, b)model_mult(a,b); % a and b are inslide the Rect()

%%  Set data path to logSNR data that we fit a model to.
data_dir_logsnr = [data_path 'included_datasets/' catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'];
data_dir_bp = [data_path 'included_datasets/' catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'];

%% Get intermodulation prominent bipolar channels. 
% Compure mean and std across trials for this channel data. 
[mean_logsnr_1, std_logsnr_1, freq_data_1] = mean_std_data(data_dir_logsnr, catname, bp_chan_1);

%% Plot a time course, log power, and logsnr for some parameters.
%% a x Rect(X) x Rect(Y)
% time course
a = 1;
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
signal_1 = a*Rect_comp1 .* Rect_comp2;

% plot time course
figure();clf
hAx = gobjects(3, 1);
hAx(1) = subplot(3, 1, 1);
plot(times, Rect_comp1);
ylabel('amplitude [-]');
title('rect(X)');
hAx(2) = subplot(3, 1, 2);
plot(times, Rect_comp2);
ylabel('amplitude [-]');
title('rect(Y)');
hAx(3) = subplot(3, 1, 3);
plot(times, signal_1);
ylabel('amplitude [-]');
title('1 x rect(X) x Rect(Y)');
xlim(hAx, [0, 0.5]);
xlabel('time [s]');
ylabel('amplitude [-]');

% power and logSNR
[logsnr_1, powers_1, freqs] = compute_logsnrs_y(signal_1, sampling_rate);
maxf = 250;
[~, maxf_ind] = find_closest(freqs, maxf);

% plot power 
figure();clf
plot(freqs(1:maxf_ind), powers_1(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logPower[dB]');

% plot logSNR
figure();clf
plot(freqs(1:maxf_ind), logsnr_1(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logSNR[-]');

%% Rect(1X) x Rect(1Y)
a = 1;
b = 1;
signal_2 = model(a, b);

% Components of models.
% f1[Hz] sine wave.
comp1 = a*sin(2*pi*f1*times);
% f2[Hz] sine wave.
comp2 = b*sin(2*pi*f2*times);
% Rectification 
Rect_comp1 = comp1;
Rect_comp1(comp1 <= 0) = 0;
Rect_comp2 = comp2;
Rect_comp2(comp2 <= 0) = 0;

% plot time course
figure();clf
hAx = gobjects(3,1);
hAx(1) = subplot(3, 1, 1);
plot(times, Rect_comp1);
ylabel('amplitude [-]');
title('rect(1X)');
hAx(2) = subplot(3, 1, 2);
plot(times, Rect_comp2);
ylabel('amplitude [-]');
title('rect(1Y)');
hAx(3) = subplot(3, 1, 3);
plot(times, signal_2);
ylabel('amplitude [-]');
title('rect(1X) x Rect(1Y)');
xlim(hAx, [0, 0.5]);
xlabel('time [s]');
ylabel('amplitude [-]');

% power and logSNR
[logsnr_2, powers_2, freqs] = compute_logsnrs_y(signal_2, sampling_rate);
maxf = 250;
[~, maxf_ind] = find_closest(freqs, maxf);

% plot power
figure();clf
plot(freqs(1:maxf_ind), powers_2(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logPower[dB]');

% plot logSNR
figure();clf
plot(freqs(1:maxf_ind), logsnr_2(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logSNR[-]');

%% Rect(1X) x Rect(0.01Y)
a = 1;
b = 0.01;

signal_3 = model(a, b);

% Components of models.
% f1[Hz] sine wave.
comp1 = a*sin(2*pi*f1*times);
% f2[Hz] sine wave.
comp2 = b*sin(2*pi*f2*times);
% Rectification 
Rect_comp1 = comp1;
Rect_comp1(comp1 <= 0) = 0;
Rect_comp2 = comp2;
Rect_comp2(comp2 <= 0) = 0;

% plot time course 
figure();clf
hAx = gobjects(3,1);
hAx(1) = subplot(3, 1, 1);
plot(times, Rect_comp1);
ylabel('amplitude [-]');
title('rect(1X)');
hAx(2) = subplot(3, 1, 2);
plot(times, Rect_comp2);
ylabel('amplitude [-]');
title('rect(0.01Y)');
hAx(3) = subplot(3, 1, 3);
plot(times, signal_3);
ylabel('amplitude [-]');
title('rect(1X) x Rect(0.01Y)');
xlim(hAx, [0, 0.5]);
xlabel('time [s]');
ylabel('amplitude [-]');

% power and logSNR
[logsnr_3, powers_3, freqs] = compute_logsnrs_y(signal_3, sampling_rate);
maxf = 250;
[~, maxf_ind] = find_closest(freqs, maxf);

% plot power
figure();clf
plot(freqs(1:maxf_ind), powers_3(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logPower[dB]');

% plot logSNR
figure();clf
plot(freqs(1:maxf_ind), logsnr_3(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logSNR[-]');


%% Rect(1X) x Rect(100Y)
a = 1;
b = 100;
signal_4 = model(a, b);

% Components of models.
% f1[Hz] sine wave.
comp1 = a*sin(2*pi*f1*times);
% f2[Hz] sine wave.
comp2 = b*sin(2*pi*f2*times);
% Rectification 
Rect_comp1 = comp1;
Rect_comp1(comp1 <= 0) = 0;
Rect_comp2 = comp2;
Rect_comp2(comp2 <= 0) = 0;

% plot time course
figure();clf
hAx = gobjects(3,1);
hAx(1) = subplot(3, 1, 1);
plot(times, Rect_comp1);
ylabel('amplitude [-]');
title('rect(1X)');
hAx(2) = subplot(3, 1, 2);
plot(times, Rect_comp2);
ylabel('amplitude [-]');
title('rect(100Y)');
hAx(3) = subplot(3, 1, 3);
plot(times, signal_4);
ylabel('amplitude [-]');
title('rect(1X) x Rect(100Y)');
xlim(hAx, [0, 0.5]);
xlabel('time [s]');
ylabel('amplitude [-]');

% power and logSNR
[logsnr_4, powers_4, freqs] = compute_logsnrs_y(signal_4, sampling_rate);
maxf = 250;
[~, maxf_ind] = find_closest(freqs, maxf);

% plot power
figure();clf
plot(freqs(1:maxf_ind), powers_4(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logPower[dB]');

% plot logsNR
figure();clf
plot(freqs(1:maxf_ind), logsnr_4(1:maxf_ind));
xlabel('freq [Hz]');
ylabel('logSNR[-]');


%% Functions
% Compute mean across trials
function [mean_data, std_data, freq_data] = mean_std_data(data_dir, catname, chan)
% Load S1 data only max condition. Then compute mean and std across trials.
%   Input:  data_dir = directory name
%           fig1_catname = cat name 
%           chan = channels we plot
%   Output: mean_data = mean across trials
%           std_data = std across trials
%           freq_data = frequencies (0-250Hz)

    % want the max stimulation condition
    fnames = dir(fullfile(data_dir, [catname, '*S1*.mat'])); 
    nConds = length(fnames);
    % load the last one (max stimulation condition)
    load(fullfile(data_dir, fnames(nConds).name))

    % Extract data from the channels (dim = frequencies x trials)
    trial_data = squeeze(data.trial(chan,:,:)); 

    % Get index of 250Hz
    [~,fend] = find_closest(data.freq{1}, 250);

    % Take mean of log SNR across trials.
    mean_data = mean(trial_data(1:fend, :), 2);
    std_data = std(trial_data(1:fend, :), 0, 2); % unbiased
    % Correct dim
    mean_data = mean_data';
    std_data = std_data';
    % Frequency data 
    freq_data = data.freq{1}(1:fend);
    freq_data = freq_data';
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

end


function model = model_plus(a, b)
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
    comp1 = a*sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = b*sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = Rect_comp1 + Rect_comp2;
end

function model = model_mult(a, b)
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
    comp1 = a*sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = b*sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = Rect_comp1 .* Rect_comp2;
end

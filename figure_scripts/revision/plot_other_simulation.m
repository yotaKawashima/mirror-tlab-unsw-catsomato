%% Model for explaining responses.
% Here, we show some models for the rebuttal letter.
% Note that square means X^2 not square wave./

%% Set variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

img_fmt = '-dpdf';

% frequency
f1 = 23;
f2 = 200;

% Sampling rate [Hz]
sampling_rate = 10000; 

% time 0s to 2s.
time_s = 0; % Start [s]
time_e = 2; % End [s]
times = time_s:1/sampling_rate:time_e;


%% Model
% Components of models.
% f1[Hz] sine wave.
comp1 = sin(2*pi*f1*times);

% f2[Hz] sine wave.
comp2 = sin(2*pi*f2*times);

% Half Wave Rectification 
%Rect_comp1 = comp1;
%Rect_comp1(comp1 <= 0) = 0;
%Rect_comp2 = comp2;
%Rect_comp2(comp2 <= 0) = 0;

% Square (x^2)
Sq_comp1 = comp1.^2;
Sq_comp2 = comp2.^2;

% Full Wave rectification |X|
FWRect_comp1 = abs(comp1);
FWRect_comp2 = abs(comp2);
%% Plot 
num_signals =7;

% Preallocate subplots.
hAx1 = gobjects(num_signals,1);    
hAx2 = gobjects(num_signals,1);

for i = 1:num_signals
    switch i 
        case 1
            % A single input frequency processed linearly.
            % X
            p_y = comp1;            
        %case 2
            % Non-linear processing of a single input. Negative values 
            % replaced with 0.
            % Rec(X)
        %    p_y = Rect_comp1;
        %case 3
            % The two inputs are first half-rectified and then summed. 
            % Rec(X) + Rec(Y)
        %    p_y = Rect_comp1 + Rect_comp2;  
        %case 4 
            % The two inputs first interact via a multiplication term and 
            % are then half rectified together. 
            % Rec(XY)
        %    p_y = comp1 .* comp2;
        %    p_y(p_y <= 0) = 0;
        %case 5
            % The two inputs are first half-rectified and then interact 
            % via a multiplication term.
            % Rec(X)*Rec(Y)
        %    p_y = Rect_comp1 .* Rect_comp2;                        
        case 2%6 
            % Square 
            % Sq(X)
            p_y = Sq_comp1;
        case 3%7 
            % Sum of two square.
            % Sq(X)+Sq(Y)
            p_y = Sq_comp1 + Sq_comp2;
        case 4%8 
            % Multyple of two square.
            % Sq(X)*Sq(Y) (=Sq(xy))
            p_y = Sq_comp1 .* Sq_comp2;
        case 5%9 
            % Full wave rectification
            p_y = FWRect_comp1;
        case 6%10 
            % Sum of two full wave rectification
            p_y = FWRect_comp1 + FWRect_comp2;
        case 7%11
            % Multyple of full wave rectification
            %|X|*|Y| (=|X*Y|)
            p_y = FWRect_comp1 .* FWRect_comp2; 
    end
    
    fig1 = figure(1);
    hAx1(i) = subplot(num_signals, 1, i);%subtightplot(9,1,i);
    p = plot(times, p_y, 'k', 'linewidth', 0.3);
    
    % Power spectrum and logSNR
    [logsnrs, powers, freqs]= compute_logsnrs_y(p_y, sampling_rate);
    
    fig2 = figure(2);
    hAx2(i) = subplot(num_signals, 1, i);%subtightplot(9,1,i);
    pp = plot(freqs, logsnrs,  'k', 'linewidth', 0.4);

    foi_f1_and_harm = 23:23:250;
    foi_f2 = 200;
    foi_inter = sort([177:-23:0 200+23:23:250]);
    foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);
    [~, ind_foi] = parfind_closest(freqs, foi);

    hold on
    % plot vertical lines at foi
    %   grab axis variables    
    ylims = get(gca, 'YLim');
    ylims = [-10, ylims(2)*0.05];
    v1 = numel(foi_f1_and_harm);
    v2 = numel(foi_inter);
    v = gobjects(numel(foi), 1);
    %   plot vertical lines
    % 23 + harmonics
    v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], ylims, ...
        'Color', [255, 47, 64]/256); 
    % 200
    v(v1+1) = plot([200, 200], ylims, ...
         'Color', [0, 205, 0]/256); 
    % intermodulation
    v(v1+2:v1+v2+1) = plot([foi_inter', foi_inter'], ylims, ...
        'Color', [0, 82, 255]/256); 
    %   format vertical lines
    set(v, 'LineWidth', 1);
    hold off
end

% Figure 1 setting
set(hAx1, 'ylim', [-1.5 1.5]) 
%set(hAx1(3), 'ylim', [-1, 2.5]);
%set(hAx1(7), 'ylim', [-1, 2.5]);
set(hAx1(3), 'ylim', [-1, 2.5]);
%set(hAx1(10), 'ylim', [-1, 2.5]);
set(hAx1(4), 'ylim', [-1, 2.5]);
%set(hAx1,'Ytick',[0])
%set(hAx1,'YTickLabel',[0]);
set(hAx1, 'xlim', [0 0.25])
set(hAx1(1:end-1), 'XTicklabel',[]);
set(hAx1,'Xtick',[0 0.125 0.25])
set(hAx1,'XtickLabel',[0 0.125 0.25])
set(hAx1(1:end-1), 'XTicklabel',[]);
xlabel(hAx1(end), 'time [s]');
%ylabel(hAx1(6), 'Amplitude [-]');
ylabel(hAx1(4), 'Amplitude [-]');
set(hAx1, 'FontSize', 8);
set(hAx2, 'FontName', 'Arial');
set(fig1,'renderer','Painters');
f=fig1;
f.Units = 'centimeters';
f.Position = [10, 10, 8, 8.5];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(fig1, img_fmt, 'FigRev2_othermodels_a');
elseif img_fmt == "-dtiff"
    print(fig1, img_fmt, 'FigRev2_othermodels_a', '-r300');
end

% Figure 2 setting
set(hAx2, 'ylim', [-10 40])
set(hAx2, 'xlim', [0 250])
set(hAx2,'Ytick',[0 30])
%set(hAx2,'YTickLabel',[0]);
set(hAx2(1:end-1), 'XTicklabel',[]);
xlabel(hAx2(end), 'frequency [Hz]');
ylabel(hAx2(4), 'logSNR [dB]');
set(hAx2, 'FontSize', 8);
set(hAx2, 'FontName', 'Arial');
set(fig2,'renderer','Painters')
f=fig2;
f.Units = 'centimeters';
f.Position = [10, 10, 8, 8.5];
if img_fmt == "-depsc" || img_fmt == "-dpdf"   
    print(fig2, img_fmt, 'FigRev2_othermodels_b');
elseif img_fmt == "-dtiff"
    print(fig2, img_fmt, 'FigRev2_othermodels_b', '-r300');
end

%% Function
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
    freqs_i = freqs(1:maxf_ind);
    
    % first frequency index at hbw*nS 
    % Some data points in the beginning and the end should be removed
    % because neighboring does not exist for them.
    % We can do zero-padding (to signal not to mask) to solve this. 
    % first frequency index at 3*nSample. 
    freqs = freqs_i(nSamp_hbw*nSamp:(end-nSamp_hbw*nSamp+1));
    nReps = length(freqs);

    logsnrs = zeros(1, nReps);

    for index = 1:nReps
        % zero pad mask
        maskm = [zeros(1, index-1) mask zeros(1,maxf_ind-mLen-index+1)];
        logsnrs(1, index) = dot(maskm, powers(1:maxf_ind));

    end
end



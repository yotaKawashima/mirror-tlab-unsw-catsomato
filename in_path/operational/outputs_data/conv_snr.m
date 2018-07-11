% data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/C20110808_R03/';
% data_name = 'C20110808_R03_TStim_S1_F023A159_F200A016_epoched_rsampsl_biprref_evkresp_cmtst3-5.mat';

%load([data_dir, data_name])

function [data_out, filt] = conv_snr(data, params)

% Uses convolution to compute the SNR.
%
% data_out = conv_snr(data, params) uses params to determine the kernal
% with which data.trial should be convolved in order to obtain the snr
% spectrum. 
%
% data_out = conv_snr(data) uses default values for params.
%
% [data_out, filt] = conv_snr(...) also returns the filter used, filt. 
%
% For dtails of params see script for default values and explanation.

if nargin < 2 || isempty(params)
    % the inner section is hardcoded to be the half bandwidth of the power
    % calculation. see TW below.
    
    params.fmax = 300; % maximum frequency to be calculated to
    
    params.kernalhbw = 3.2; 
    % the half-length of the kernal. this includes the central part as well
    % as any zeros.
    
    params.kernalzerolen = 1/3;
    % Any zeros to be padded between the inner high and the outer low. Set
    % to 0 for no padding. 
end


% copy everything
data_out = data;

% determine trial length, T
T = data.custom.time(end) - data.custom.time(1);

% final frequency = 300;
f_ind = find_closest_ind(data.freq{1}, params.fmax);

% pull TW
TW = data.custom.params.tapers(1);
W = TW/T; % half bandwidth

% so we need to start that far out and then go out to f + F
F = params.kernalhbw;

% determine the indices that we need. 
f_res = mean(diff(data.freq{1}));
inner = round(W/f_res);
outer = round(F/f_res);
middle = round(W*(1+params.kernalzerolen)/f_res);

% set up the convolution filter
f_in = ones(1,2*inner+1);
f_0 = zeros(1, middle-inner);
f_out_tmp = -ones(1,outer-middle);
f_out_total = sum(f_in);
f_out= f_out_tmp*f_out_total/2/numel(f_out_tmp);
filt = [f_out, f_0, f_in, f_0, f_out];


data_in = data.trial(:, 1:f_ind, :);
[nCh, nSamp, nTr] = size(data_in);

trialtmp = zeros(nCh, nSamp+length(filt)-1, nTr);

for cCh = 1:nCh
    for k = 1:nTr

        trialtmp(cCh, :, k) = conv(data_in(cCh, :, k), filt);
        
    end
end

% cleaning, like Levi from Attack on Titan
data_out.trial = trialtmp;
data_out.freq{1} = data.freq{1}(1:f_ind);
data_out.custom.convkernal = filt;

if nargout < 2
    clear filt
end
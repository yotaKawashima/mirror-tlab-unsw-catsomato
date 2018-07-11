function data_out = snrtosurrounds(data, maxf, innerdiff, outerdiff)

% snrtosurrounds: calculates the SNR with respect to neighbours. 
%   data_out = snrtosurrounds(data, maxf, innerdiff, outerdiff) calculates
%   the SNR for the frequencies between 0 and maxf. innerdiff and outerdiff
%   define the neighbourhood region (see below).
% 
%   data_out = snrtosurrounds(data) uses the default values of maxf=250,
%   innerdiff=1, outerdiff=3.
% 
%   The neighbourhood region is illustrated in the following diagram.
%   F=frequency of interest
%   A=F +- outerdiff
%   B=F +- innerdiff
%   C=F +- halfbandwidth (0.5Hz)
%   
%              C_F_C
%    __A   B___|   |___B   A__
%      |___|           |___|


% Written by Rannee Li, Jun 2017


if nargin < 2 % set defaults
    innerdiff = 1;
    outerdiff = 3;
    maxf = 250;
end

data_out = data;
data_out.trial = [];

% rename
nT = data.custom.ntrials;
nSig = data.custom.nsignals;

% create the mask
% first, find the number of indicies I need to be 1 Hz
% assume they're all evenly spaced. 
spacing = mean(diff(data.freq{1}));
nSamp = round(1/spacing);

subFl = outerdiff-innerdiff;
subFh = innerdiff;

% mask
mLen = 2*(2*floor(nSamp/2)+nSamp*2)+1;
mask = zeros(1, mLen);

% neighbour mask
lomask = 1/(4*nSamp);
mask(1:nSamp*2) = -lomask;
mask(end-nSamp*2+1:end) = -lomask;

% central mask
himask = 1/(floor(nSamp/2)*2 + 1);
m = ceil(mLen/2);
mask(m-3:m+3) = himask;

% expand and pad for vectorisation
mask = repmat(mask, [nSig,1,nT]);


% reduce to only frequencies of interest
[~, k] = find_closest(data.freq{1}, maxf);
freq_i = data.freq{1}(1:k);

% first frequency index at 3*nS
data_out.freq{1} = freq_i(3*nSamp:(end-3*nSamp+1));
nReps = length(data_out.freq{1});

% disp(nReps)

snrs = zeros(nSig, nReps, nT);

%tic;

for m = 1:nReps
    % pad mask
    maskm = [zeros(nSig, m-1, nT) mask zeros(nSig, k-mLen-m+1, nT)];
    
    snrs(:, m, :) = sum(maskm.*data.trial(:, 1:k, :), 2);
    
%     if mod(m,10)==0
%         disp([num2str(m, '%- 5i') '  ' num2str(toc, '%.2f')])
%     end
    
end


data_out.trial = snrs;
clear snrs

% modify file name
data_out.custom.filename = [data.custom.filename '_snrsurr'];
data_out.custom.snrbounds = [innerdiff outerdiff];
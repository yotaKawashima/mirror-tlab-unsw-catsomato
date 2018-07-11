data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/C20110808_R03/';
data_name = 'C20110808_R03_TStim_S1_F023A159_F200A016_epoched_rsampsl_biprref_evkresp_cmtst3-5.mat';

%load([data_dir, data_name])

% determine trial length, T
T = data.custom.time(end) - data.custom.time(1);

% final frequency = 300;
fend = 300;
f_ind = find_closest_ind(data.freq{1}, fend);

% pull TW
TW = data.custom.params.tapers(1);
W = TW/T; % half bandwidth

% so we need to start that far out and then go out to f + F
F = 3.2;

% determine the indices that we need. 
f_res = mean(diff(data.freq{1}));
inner = round(W/f_res);
outer = round(F/f_res);
middle = round(4*W/3/f_res);

% set up the convolution filter
filt_half = [-ones(1, outer-middle)/(outer-middle), zeros(1, middle-inner), ones(1, inner)/(inner)];
filt = [filt_half fliplr(filt_half)];


data_in = data.trial(256, 1:f_ind, :);
[nCh, nSamp, nTr] = size(data_in);

trialtmp = zeros(nCh, nSamp, nTr);

for cCh = 1:nCh
    for k = 1:nTr

        trialtmp(cCh, :, k) = conv(data_in(cCh, :, k), filt, 'same');



    end
end
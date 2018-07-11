% load data
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/epoched_rsampsl_biprref_evkresp_cmtspwr_adatout';

catname = 'C20110808_R03';

fname = dir(fullfile(data_dir, [catname '*S1' '*.mat'])); 

%load(fullfile(data_dir, fname(1).name))


%% P-value processing
q = 0.05;
threshval = 'default';

% To classify channels at a given frequency
plotopts(5).singfoi      = true;  % DO NOT EDIT
plotopts(5).singchan     = false; % DO NOT EDIT
plotopts(5).classify     = true;  % DO NOT EDIT
% foi (below) are the frequencies of interest. Each frequency has its own
% plot.
plotopts(5).foi = [23, 200];
% classlim (below) is the p-value for which p<classlim is significant.
plotopts(5).classlim = threshval;
plotopts(5).fdr = {q, 'pID'};

anp_figures_ph(data_dir, [catname '*S1'], plotopts(5))
options = plotopts;

rep = 5;
q = options(rep).fdr{1};
[pID,pN] = FDR(pvals,q);
eval(['options(rep).classlim = ' options(rep).fdr{2} ';']);

%%
p_lim = 0.05; %options(5).classlim;

% first classify the channels
p_thresh = zeros(size(pvals));
p_thresh(pvals<p_lim) = 1;
% convert binary to decimal - int/(:, :, 3) is LSB
p_class = p_thresh;
p_class(:, :, 2) = p_class(:, :, 2)*2;
p_class(:, :, 1) = p_class(:, :, 1)*4;
p_class = sum(p_class, 3);

% pie chart by p-value
f23_i = find_closest_ind(metavars.freq{1}, 23);
p_resp23 = p_class(:, f23_i);
p_resp23 = p_resp23(p_resp23>0);

[N,EDGES] = histcounts(p_resp23);
pie(N)
% this is the script version of the interactions plot. 

%% shared variables
% beginning of the filename that contains the data to be analysed
filename_out = 'C20110808_R01';

% the directory where the auxiliary functions can be found
% generally does not beed to be edited.
% func_dir = '/home/rannee/Documents/MATLAB/aux_files/';
func_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/';

% the directory where the data can be found
% generally does not beed to be edited.
% data_dir = '/home/rannee/Documents/MATLAB/data/';
data_dir = ['/media/phoebeyou/My Passport/Spencers_Cat_Data/' filename_out '/'];

% the directory where the figures should be put
% generally does not beed to be edited.
% fig_dir = '/home/rannee/Documents/MATLAB/figures/';
fig_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';

% type of z-scoring
% 1 = zscore, 2 = new_zscore, 0 = none
z_type = 0; 

% did you use chronux to find power spectrum?
chronspec = true;

% areas used in analysis
A = {'_S1', '_S2'};

% frequencies of interest for anova analysis
% FOI = 31; %[20, 23, 46, 69, 154, 177, 200];

% the p-value, threshval, for which p < threshval is significant
% threshval = 4.7922e-05;

% the value associated with an error to do with moving the files
sys = 1; % if on mac, this should be 64


% and onto other orders of business
% add the functions to path - do not edit
addpath(genpath(func_dir))

% assign the z-score directory - do not edit
if z_type == 1
    z_dir = '_zscored';
elseif z_type == 2
    z_dir = '_newzscr';
else
    z_dir = [];
end
if chronspec
    p_dir = 'cmtspwr';
else
    p_dir = 'pwrspec';
end


%% load data
dir_here = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];
%fname = dir(fullfile(dir_here, [filename_out '*.mat']));
% this gives us all of the conditions though. we have to load the first one
% to work out how many there are and what their names are

% for now, let's work with constant 20/23 and varying 200, so we'll
% hardcode that in.
fname = dir(fullfile(dir_here, [filename_out '*S1*F023A079*.mat']));

fprintf('Loading data ')
for i = 1:length(fname)
    load([dir_here fname(i).name])
    datas(i) = data;
    fprintf('.')
end

fprintf('\n')

%% plot
% I remember the effect was most distinct at 23 Hz
f1 = str2double(data.custom.filename(25:27));
f2 = str2double(data.custom.filename(34:36));

f_stim = sort([f1, f2]);

[~, fl_ind] = find_closest(data.freq{1}, f_stim(1));

plotmat = zeros(data.custom.nsignals, data.custom.subplotconfig(1), data.custom.ntrials);

for i = 1:data.custom.subplotconfig(2)
    plotmat(:, i, :) = datas(i).trial(:, fl_ind, :);
    
end

x = mean(plotmat, 3);

figure()
%x1 = squeeze(plotmat(1:100,:,1))';
x1 = x';
plot(x1)

% center this plot at zero.
basex = x1(1, :);
basex = repmat(basex, data.custom.subplotconfig(2), 1);

x2 = x1-basex;
l = find_lims(x2);
figure(1); clf
plot(x2(:, 1:50))
set(gca, 'YLim', l)
figure(2); clf
plot(x2(:, 51:100))
set(gca, 'YLim', l)






































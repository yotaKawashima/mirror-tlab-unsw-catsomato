% This m-file calls the functions that are needed in order to perform
% anova2 analysis on the data. it requires that the data is pre-processed
% into the correct format.
%
% Edit the first section to select the cat. You may need to run this
% section even if you do not intend to run all other sections.

% Last edited by Rannee Li, Jul 2015
% Change by Phoebe You, 12 Jan 2016

clear all;
close all;
clc;

dbstop if error

Job.prepare = 0;
Job.process = 1;

% Define which frequency range used for computing 10,20,30,40,50?
FreqSeg = 10;

Job.plot = 0;
q = 0.05; %0.01; % q=0.1 --- use all experiments together to get 1 pthr. 

%% pre-preprocessing - EDIT THIS SECTION
% beginning of the filename that contains the data to be analysed
filename_out = 'C20110808_Rx4_TStim';
disp(filename_out);

% the directory where the auxiliary functions can be found
% generally does not beed to be edited.
% func_dir = '/home/rannee/Documents/MATLAB/aux_files/';
func_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/';

% the directory where the data can be found
% generally does not beed to be edited.
% data_dir = '/home/rannee/Documents/MATLAB/data/';
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/raw/';

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
%% prepare the data
if Job.prepare
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];

TargFreq = 100;
for area = 1:numel(A)
filename_header = [filename_out A{area}];
%% frange = [10 5] --- [segment values, number of segment pairs] eg,get the mean of 90-110,80-120,70-130,60-140,50-150
% the segment value is 10, and the number of segment pairs is 5
data = anp_data_around100(dir_sec, filename_header, [10 5], false, TargFreq);
%% Use for further analysis about 200+n*23 response
% data = anp_data(dir_sec, filename_header, [0 300], false);
end

mverr = system(['mv ' filename_out '*.mat ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/']);

if mverr == sys
    warning('Destination directory not found. Making directory.')
    system(['mkdir ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/']);
    system(['mv ' filename_out '* ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/']);
end
end
%% do the analysis
if Job.process
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/'];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

% anp_analysis_around100(dir_sec, filename_header, FreqSeg);
anp_analysis(dir_sec, filename_header);

end

mverr = system(['mv ' filename_out '*.mat ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatout/']);

if mverr == sys
    warning('Destination directory not found. Making directory.')
    system(['mkdir ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatout/']);
    system(['mv ' filename_out '* ' data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatout/']);
end
end
%% draw figures
if Job.plot
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatout/'];

% chan{1} = [178, 146];
% chan{2} = [26, 45];
chan{1} = [143, 257];
chan{2} = [172, 165];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

% threshval = 25;
threshval = 'default';

% To plot p-values for all frequencies and all channels
plotopts(1).singfoi      = false; % DO NOT EDIT
plotopts(1).singchan     = false; % DO NOT EDIT
plotopts(1).classify     = false; % DO NOT EDIT
% clim (below) is the max value of the colorbar. set to [] if you want to
% use default.
plotopts(1).clim = threshval;    


% To classify all channels and all frequencies                                
plotopts(2).singfoi      = false; % DO NOT EDIT
plotopts(2).singchan     = false; % DO NOT EDIT
plotopts(2).classify     = true;  % DO NOT EDIT
% classlim (below) is the p-value for which p<classlim is significant.
% plotopts(2).classlim = threshval;
plotopts(2).classlim = 'fdr';
plotopts(2).fdr = {q, 'pID'};
% plotopts(2).fdr = {q};

% To plot p-values for all frequencies for a single channel
plotopts(3).singfoi      = false; % DO NOT EDIT
plotopts(3).singchan     = true;  % DO NOT EDIT
plotopts(3).classify     = false; % DO NOT EDIT
% clim (below) holds the values of horizontal dotted lines (p-values of
% interest). If you don't want them, set to [].
plotopts(3).clim = [0.01, 0.001, threshval];
% foi (below) holds the values of the vertical dotted lines (frequencies of
% interest). If you don't want them, set to [].
plotopts(3).foi = [23, 46, 200];
% plotopts(3).foi = [];
% chan (below) holds the channels. Each channel appears on a separate plot.
plotopts(3).chan = chan{area};

% To plot p-values for all channels at a given frequency
plotopts(4).singfoi      = true;  % DO NOT EDIT
plotopts(4).singchan     = false; % DO NOT EDIT
plotopts(4).classify     = false; % DO NOT EDIT
% clim (below) is the max value of the colorbar. set to [] if you want to
% use default.
plotopts(4).clim = threshval;
% foi (below) are the frequencies of interest. Each frequency has its own
% plot.
plotopts(4).foi = [23, 200];

% To classify channels at a given frequency
plotopts(5).singfoi      = true;  % DO NOT EDIT
plotopts(5).singchan     = false; % DO NOT EDIT
plotopts(5).classify     = true;  % DO NOT EDIT
% foi (below) are the frequencies of interest. Each frequency has its own
% plot.
plotopts(5).foi = [23, 200];
% classlim (below) is the p-value for which p<classlim is significant.
plotopts(5).classlim = threshval;

anp_figures_ph(dir_sec, filename_header, plotopts(2))

end

mverr = system(['mv ' filename_out '*.png ' fig_dir 'anovafigures/']);

if mverr == sys
    warning('Destination directory not found. Making directory.')
    system(['mkdir ' fig_dir 'anovafigures/']);
    system(['mv ' filename_out '*.png ' fig_dir 'anovafigures/']);
end
end
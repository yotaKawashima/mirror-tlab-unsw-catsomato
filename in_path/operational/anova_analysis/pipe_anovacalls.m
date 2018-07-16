% This m-file calls the functions that are needed in order to perform
% anova2 analysis on the data. it requires that the data is pre-processed
% into the correct format.
%
% Edit the first section to select the cat. You may need to run this
% section even if you do not intend to run all other sections.

%% pre-preprocessing - EDIT THIS SECTION
% beginning of the filename that contains the data to be analysed
filename_out = 'C20110808_R03';

% add directories to path
run(fullfile(mfilename('fullpath'), '../../../../figure_scripts/path_setup.m'))

% do you wish to use chronux to calculate power spectrum?
chronspec = true;

% the directory where the data can be found
% generally does not need to be edited.
data_dir = fullfile(data_path, 'included_datasets', filename_out, '/');

% type of z-scoring
% 1 = zscore, 2 = new_zscore, 0 = none
z_type = 0; 

% areas used in analysis
A = {'_S1', '_S2'};

% and onto other orders of business

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
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

data = anp_data(dir_sec, filename_header, [0 220], false);

end

afpc_filemover(filename_out, dir_sec, 'adatain')
%% do the analysis
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/'];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

anp_analysis(dir_sec, filename_header);

end

afpc_filemover(filename_out, dir_sec, 'adatout')

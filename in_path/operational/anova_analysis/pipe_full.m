% This m-file calls the functions that are needed in order to process the
% raw data into power spectrum data that can be fed into further analysis. 
%
% Edit the first section to select the cat. You may need to run this
% section even if you do not intend to run all other sections.

%% options and selections - EDIT THIS SECTION
% filename of the RAW data that is to be used. 
filename_in = 'Cat20110808_UtahRight-3_LFPStimSegs';

% filename to be saved. 
filename_out = 'C20110808_R03';

% add directories to path
run(fullfile(mfilename('fullpath'), '../../../../figure_scripts/path_setup.m'))

% do you wish to use chronux to calculate power spectrum?
chronspec = true;

% the directory where the data can be found
% generally does not beed to be edited.
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

%% preprocessing - 1: epoch the data.
% call epoch_and_separate
load(fullfile(data_dir, '../../raw/', [filename_in '.mat']))

output_dest = [data_dir 'epoched/'];

system(['mkdir ' data_dir]);
system(['mkdir ' output_dest]);

epoch_and_separate(filename_in, tres, timestamp, lfp1_data, lfp2_data, ...
    stimFreqY, stimFreqX, stimMeshY, stimMeshX, stimTime, output_dest, filename_out)

clear tres timestamp lfp1_data lfp2_data stimFreqY stimFreqX stimMeshY stimMeshX
%% preprocessing - 2: resample
% call resampler
dir_sec = [data_dir 'epoched/'];

loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));

stimDur = stimTime.rampup + stimTime.presine + stimTime.sine + stimTime.postsine + stimTime.rampdown;
target = -0.5:0.0001:stimDur+0.5; % this is the target in terms of the original timestamp
target = target - stimTime.rampup - stimTime.presine;

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    data = resample_chooser(data, 'spline', target);
    
    save(data.custom.filename, 'data')
    
    clear data
end

afpc_filemover(filename_out, dir_sec, 'rsampsl')

%% preprocessing - 3: bipolar rereferencing
% call bipolreref
dir_sec = [data_dir 'epoched_rsampsl/'];

loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    data = bipolreref(data);
    
    save(data.custom.filename, 'data')
    
    clear data
end

afpc_filemover(filename_out, dir_sec, 'biprref')
%% preprocessing - 4: extract evoked response0..
% call extractstableresponse - for timepoints

dir_sec = [data_dir 'epoched_rsampsl_biprref/'];

loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    data = extractstableresponse(data, [0.5 2.5]);
    
    
    save(data.custom.filename, 'data')
end

afpc_filemover(filename_out, dir_sec, 'evkresp')

%% preprocessing - 5: zscore
if z_type
    dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp/'];
    
    loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));
    
    for i = 1:numel(loadname)
        load(fullfile(dir_sec, loadname(i).name))
        
        if z_type == 1
            data = zscore(data);
        elseif z_type == 2
            data = new_zscore(data);
        end
        
        save(data.custom.filename, 'data')
    end
    
    afpc_filemover(filename_out, dir_sec, z_dir(2:end))

end
%% preprocessing - 6: compute power spectrum
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '/'];

loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));

if chronspec
    addpath(genpath(chron_dir))
end

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    if chronspec
        data = chronux_pwrspec(data);
    else
        data = pwrspec(data); 
    end
    
    save(data.custom.filename, 'data')
end

afpc_filemover(filename_out, dir_sec, p_dir)

%% anova: prepare the data
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

data = anp_data(dir_sec, filename_header, [0 220], false);

end

afpc_filemover(filename_out, dir_sec, 'adatain')
%% anova: do the analysis
dir_sec = [data_dir 'epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain/'];

for area = 1:numel(A)
filename_header = [filename_out A{area}];

anp_analysis(dir_sec, filename_header);

end

afpc_filemover(filename_out, dir_sec, 'adatout')

%% SNR: Find SNR
data_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];

call_snrtosurrounds(data_path, data_sec, filename_out)

%% anova: grab the fstats only
% this calculates it for all known cats, skipping those that have already
% been done.
data_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain_adatout/'];

extract_fstat(data_path, data_sec)

%% HGP: calculate the evoked power
data_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '/'];

call_evokedpower(data_path, data_sec, filename_out)

%% HGP: find the total power in the HG band
data_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_evkdpwr/'];

call_highgammapower(data_path, data_sec, filename_out)

%% HGP: anova
dir_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_evkdpwr_hgpcomp/'];

call_hgpanova(data_path, dir_sec, filename_out)

%% HGP: anova - grab the fstats only
% this calculates it for all known cats, skipping those that have already
% been done.
data_sec = ['epoched_rsampsl_biprref_evkresp' z_dir '_' p_dir '_adatain_adatout/'];

extract_fstat(data_path, data_sec)



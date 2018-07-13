% This m-file calls the functions that are needed in order to process the
% raw data into power spectrum data that can be fed into further analysis. 
%
% Edit the first section to select the cat. You may need to run this
% section even if you do not intend to run all other sections.

%% options and selections - EDIT THIS SECTION
% filename of the RAW data that is to be used. 
filename_in = 'Cat20090505_Utah-21_LFPStimSegs';

% filename to be saved. 
filename_out = 'C20110511_R02';

run('../../../figure_scripts/path_setup.m')

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

% the value associated with an error to do with moving the files
sys = 1; % if on mac, this should be 64


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
load([data_dir '../raw/' filename_in '.mat'])

output_dest = [data_dir 'epoched/'];

system(['mkdir ' output_dest]);

epoch_and_separate(filename_in, tres, timestamp, lfp1_data, lfp2_data, stimFreqY, stimFreqX, stimMeshY, stimMeshX, stimTime, output_dest, filename_out)

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

afpc_filemover(filename_out, dir_sec, 'rsampsl', sys)

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

afpc_filemover(filename_out, dir_sec, 'biprref', sys)
%% preprocessing - 4: extract evoked response0..
% call extractstableresponse - for timepoints

dir_sec = [data_dir 'epoched_rsampsl_biprref/'];

loadname = dir(fullfile(dir_sec, [filename_out, '*.mat']));

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    data = extractstableresponse(data, [0.5 2.5]);
    
    
    save(data.custom.filename, 'data')
end

afpc_filemover(filename_out, dir_sec, 'evkresp', sys)

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
    
    afpc_filemover(filename_out, dir_sec, z_dir(2:end), sys)

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

afpc_filemover(filename_out, dir_sec, p_dir, sys)

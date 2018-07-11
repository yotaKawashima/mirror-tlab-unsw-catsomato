%% Timing test 1: fullfile and if, or for and fullfile

% test 1
clear

tic;

recordings = {'C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'};

data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';

area = 'S1';

for k = 1:numel(recordings)
    fname = dir(fullfile(data_dir, [recordings{k} '_' area '*.mat']));
    disp(fname.name)
end



time1 = toc;
disp(time1)


% test 2
clear

tic;

recordings = {'C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'};

data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';

area = 'S1';

fname = dir(fullfile(data_dir, recordings));


for k = 1:numel(fname)
    if strfind(fname(k).name, area)
        disp(fname(k).name)
    end
end

time1 = toc;
disp(time1)

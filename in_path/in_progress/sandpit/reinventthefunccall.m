%% the files got rearranged and my old ways don't work no more D:

catname = 'C20110808_R01';

addpath(genpath('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files'))
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';

% dir_sec = [data_dir catname '/epoched_rsampsl_biprref_evkresp/']; %_dwnsmpl

chron_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/thirdparty_toolboxes/chronux';


%% cmtspwr
dir_sec = [data_dir catname '/epoched_rsampsl_biprref_evkresp/'];

loadname = dir(fullfile(dir_sec, [catname, '*.mat']));


addpath(genpath(chron_dir))


for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))

    data = chronux_pwrspec(data);
    
    save(data.custom.filename, 'data')
end

afpc_filemover(catname, dir_sec, 'cmtspwr', 1)

rmpath(genpath(chron_dir))
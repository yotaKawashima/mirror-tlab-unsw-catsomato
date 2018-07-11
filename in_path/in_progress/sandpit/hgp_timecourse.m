% HGP timecourse

% load data
recordings = {'C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'};

area = {'S1', 'S2'};

data_dir_lvl1 = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';

data_dir_lvl3 = '/epoched_rsampsl_biprref/';

data_dir_lvl2 = 'raw/';

nRec = numel(recordings);

addpath(genpath('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/thirdparty_toolboxes/chronux'))

for cRec = 5%1:nRec
    
    % realign time
    %fname_og = dir(fullfile(data_dir_lvl1, data_dir_lvl2, [recordings{cRec} '*' area{cArea} '*.mat']));
    %load(fullfile(data_dir,fname.name))
    
    for cArea = 1:2
        
        % load data
        fname = dir(fullfile(data_dir_lvl1, recordings{cRec}(1:13), data_dir_lvl3, [recordings{cRec} '*' area{cArea} '*.mat']));
        
        for cCond = 1:numel(fname)
        
            load(fullfile(data_dir_lvl1, recordings{cRec}(1:13), data_dir_lvl3, fname(cCond).name))
            
            % how to break into sections? Maybe use the 4th dimension? we'd
            % need to reallocate it to work with the power spectrum anyway
            % or maybe the best solution is to actually do the power
            % spectrum and then save it all together. 
            
            % update legacy data, just in case
            data.time = data.time(1);
            
            % need to do this otherwise we overwrite
            data_tmp = data;
            
            % take this out separately and remove the unipolar channels.
            nunipol = prod(data.custom.spatialconfig);
            trialtmp = data.trial((nunipol+1):end, :, :);
            data_tmp.custom.nsignals = size(trialtmp, 1);
            
            
            for k = 1:floor(data.custom.nsamples/5000)
                fprintf('%i\t', k)
                
                data_tmp.trial = trialtmp(:, ((k-1)*5000+1):(k*5000), :);
                data = chronux_pwrspec(data_tmp);
                
                save([data_dir_lvl1, 'data/HGP_data/', fname(cCond).name(1:15), ...
                    num2str(k, '%03i'), 'T', fname(cCond).name(20:end-4), '_hgpower'], 'data')
                
                
            end 
            
            clear data data_tmp
            
        end

        
        % do ANOVA for each slice of time
        % load
        
        % data = anp_data_v2(data_dir, filename_header, [], 0);
        
        
        % visualise
        
        
    end
end

rmpath(genpath('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/thirdparty_toolboxes/chronux'))

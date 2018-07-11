catname = 'C20110808_R03';

savename_head = [catname '_TStim_'];
savename_foot = 'epoched_rsampsl_biprref_evkresp_cmtspwr';

dir_sec = ['/media/phoebeyou/My Passport/Spencers_Cat_Data/data/' savename_foot];

areas = {'S1'; 'S2'};

for a = 1:2
    fname = dir(fullfile(dir_sec, [catname, '*', areas{a}, '*.mat'])); 
    
    for f = 1:numel(fname)
        load(fullfile(dir_sec, fname(f).name))
        alltrial(:, :, :, f) = data.trial; 
    end
    
    save([savename_head areas{a} '_' savename_foot '_alltrial'], 'alltrial')
    
end % end for on a






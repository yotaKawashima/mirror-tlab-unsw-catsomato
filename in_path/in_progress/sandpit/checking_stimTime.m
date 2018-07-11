load('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/allcat_meta.mat')

nRec = numel(datas);

stim_dat = zeros(nRec, 6);

for k = 1:numel(nRec)
    
    stimtmp = datas(k).stimTime;
    
    stim_dat(k, 2) = stimtmp.rampup;
    stim_dat(k, 3) = stimtmp.presine;
    stim_dat(k, 4) = sine;
    stim_dat(k, 5) = 
    
    
    
    
end
filename = {'C20110511_R02';'C20110510_R05';...
    'C20110510_R06';'C20110808_R01';'C20110808_R03';...
    'C20110808_R04';'C20110808_R06';'C20110808_R09';'C20110808_Rx4'};
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';


nRec = numel(filename);



%%
input_type = '/epoched_rsampsl_biprref_evkresp_cmtst3-5/';
out_type = 'adatain';

for a = 1%[1, 2]
for m = 5:numel(filename)
    dir_sec = [data_dir filename{m} input_type];
    
    loadname = dir(fullfile(dir_sec, [filename{m}, '*S', num2str(a), '*.mat']));
    
    
    
    for k = 1:numel(loadname)
        %load(fullfile(dir_sec, loadname(k).name))
        
        
        data = anp_data_v2(dir_sec, [filename{m}, '*S', num2str(a)], [0 300], 2);
        
        
    end
    
    %afpc_filemover(filename{m}, dir_sec, out_type, 1)
    
end

end



%%
input_type = 'epoched_rsampsl_biprref_evkresp_cmtst3-5_adatin/';

for a = [1, 2]
for m = 1:numel(filename)
    dir_sec = [data_dir];% filename{m} input_type];
    
    %loadname = dir(fullfile(dir_sec, [filename{m}, '*S', num2str(a), '*.mat']));
    
    
    
 %   for k = 1:numel(loadname)
        %load(fullfile(dir_sec, loadname(k).name))
        
        
        anp_analysis_2(dir_sec, [filename{m}, '*S', num2str(a)])
        
        
  %  end
    
    %afpc_filemover(filename{m}, dir_sec, out_type, 1)
    
end

end

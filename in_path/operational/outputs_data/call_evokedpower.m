function call_evokedpower(data_path, datatype_in, cat_name)
% Could not find this file in Rannee's repo. 
% So, write this by myself. Maybe result in errors. 
% Written by Yota, Feb 2020.

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

% separate by S1 and S2
for a = 1:2
    
    loadname = dir(fullfile(dir_sec, [cat_name '*S' num2str(a) '*.mat']));
    
    for i = 1:length(loadname)
        
        load(fullfile(dir_sec, loadname(i).name))
        if i == 1 % Baseline condition : Amplitude of F1, F2 = (0,0)
            % Compute baseline
            baseline = calc_evoked_power(data);
            
            % Compute vibration evoked log power.
            data = calc_evoked_power(data, baseline, 0);
        else
            % Compute vibration evoked log power.
            data = calc_evoked_power(data, baseline, 0);
        end
        
        % Copy file into *_evkdpwr dir.
        save(data.custom.filename, 'data')
              
    end
    clear baseline
end

afpc_filemover(cat_name, dir_sec, 'evkdpwr')
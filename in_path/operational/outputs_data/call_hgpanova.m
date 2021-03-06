function call_hgpanova(data_path, datatype_in, cat_name)

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

% find data
loadname = dir(fullfile(dir_sec, [cat_name '*.mat']));

for a = 1:2
    
    load(fullfile(dir_sec, loadname(a).name))

    %hgp_anova(data, dir_sec)
    hgp_anova_y(data, dir_sec)
    
    % Extract only fstats 
    data_path_a = [dir_sec(1:end-1) '_adatain_adatout/'];
    extract_fstat_y(data_path_a);
    afpc_filemover(cat_name, data_path_a, 'fstonly');

end


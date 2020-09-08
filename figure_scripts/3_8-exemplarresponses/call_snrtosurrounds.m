function call_snrtosurrounds(data_path, datatype_in, cat_name)

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

loadname = dir(fullfile(dir_sec, [cat_name, '*.mat']));

for i = 1:numel(loadname)
    load(fullfile(dir_sec, loadname(i).name))
    data = snrtosurrounds(data);
    
    save(data.custom.filename, 'data')
    
    clear data
end

afpc_filemover(cat_name, dir_sec, 'snrsurr')
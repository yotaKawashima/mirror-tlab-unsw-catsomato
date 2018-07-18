function call_evokedpower(data_path, datatype_in, cat_name)

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

for a = 1:2

    loadname = dir(fullfile(dir_sec, [cat_name, '*S' num2str(a) '*.mat']));

    load(fullfile(dir_sec, loadname(1).name))
    baseline = data;

    for i = 2:numel(loadname)
        load(fullfile(dir_sec, loadname(i).name))
        data = calc_evoked_power(data, baseline, []);

        save(data.custom.filename, 'data')

        clear data
    end

    afpc_filemover(cat_name, dir_sec, 'evkdpwr')

end
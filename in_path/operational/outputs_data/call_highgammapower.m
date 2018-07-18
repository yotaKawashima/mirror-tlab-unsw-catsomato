function call_highgammapower(data_path, datatype_in, cat_name)

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

% separate by S1 and S2
for a = 1:2
    loadname = dir(fullfile(dir_sec, [cat_name '*S' num2str(a) '*.mat']));

    data = highgammapower(dir_sec, loadname, [50 150], foi, 0.5);

    save(data.custom.filename, 'data')

end

afpc_filemover(cat_name, dir_sec, 'hgpcomp')
function call_highgammapower(data_path, datatype_in, cat_name)

data_dir = fullfile(data_path, 'included_datasets', cat_name, '/');
dir_sec = fullfile(data_dir, datatype_in, '/');

% create foi list
allfoi = [23:23:150, 200:-23:50];
allfoi = sort(unique(allfoi));
foi = allfoi(allfoi<150);
foi = foi(foi>50);

% separate by S1 and S2
for a = 1:2
    loadname = dir(fullfile(dir_sec, [cat_name '*S' num2str(a) '*.mat']));
    
    % Just extract VELogP data in 50Hz<f<150Hz. (Frequencies of interest
    % are not included.)
    %data = highgammapower(dir_sec, loadname, [50 150], foi, 0.5);
    data = highgammapower_y(dir_sec, loadname, [50 150], foi, 0.5);
    
    % Save extracted data. (Mean is not taken across frequencies and trials
    % at this point.)
    save(data.custom.filename, 'data')

end

afpc_filemover(cat_name, dir_sec, 'hgpcomp')
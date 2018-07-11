function names = dirsinside(data_dir)

% finds all directories one level down


d = dir(data_dir);
isub = [d(:).isdir]; %# returns logical vector
names = {d(isub).name}';
names(ismember(names,{'.','..'})) = [];
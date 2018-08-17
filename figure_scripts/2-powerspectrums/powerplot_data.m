function data = powerplot_data(data_dir, data_type, filename, savedata)

if nargin < 3
    filename = [];
end
if nargin < 4
    savedata = false;
end

% find out if the file exists. 
if ~isempty(filename)
    if exist(filename, 2)
        load(fullfile(filename))
        return
    end
end


data_out = [];

cat_names = dirsinside(data_dir);

for c = 1:numel(cat_names)
    
    loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, '*.mat'));
    
    for k = 1:numel(loadnames)
        % load
        load(fullfile(data_dir, cat_names{c}, data_type, loadnames(k).name))
        
        % get it into the avtrcon format
        data_tmp = data.trial;
        data_tmp = mean(data.trial, 3);
        
        % smush
        data_out = [data_out;data_tmp];
    end
end

data_out(isinf(data_out)) = NaN;
data.trial = nanmean(data_out, 1);
data.custom.filename(2:13) = 'xxxxxxxx_Rxx';

if nargout < 1 || savedata
    % save this data
    save(data.custom.filename, 'data')
end

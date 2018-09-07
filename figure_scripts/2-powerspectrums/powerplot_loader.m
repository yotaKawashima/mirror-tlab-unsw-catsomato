function data = powerplot_loader(data_dir, data_type, save_dir, cond, dontsave)

% data_dir = directory directly above the one containing the directories
% for the cats

if nargin < 4
    cond = '*';
end

% find out if the file already exists
fname = dir(fullfile(save_dir, ['*' data_type '*.mat']));
if ~isempty(fname)
    load(fullfile(save_dir, fname.name))
    return
end


% find names of all cats
cat_names = dirsinside(data_dir);

data_out = [];

for c = 1:numel(cat_names)
    % check if the necessary files exist
    
    
    % find the names of all files
    loadnames = dir(fullfile(data_dir, cat_names{c}, data_type, ['*' cond '*.mat']));
    
    for k = 1:numel(loadnames)
        % load
        load(fullfile(data_dir, cat_names{c}, data_type, loadnames(k).name))
        
        % remove the third dimension
        data_tmp = [];
        for d = 1:size(data.trial, 3)
            data_tmp = [data_tmp; data.trial(:, :, d)];
        end
               
        % smush
        data_out = [data_out;data_tmp];
    end
end

data_out(isinf(data_out)) = NaN;
data.trial = nanmean(data_out, 1);
data.custom.filename(2:13) = 'xxxxxxxx_Rxx';

if nargout < 1 && ~dontsave
    % save this data
    save(data.custom.filename, 'data')
end


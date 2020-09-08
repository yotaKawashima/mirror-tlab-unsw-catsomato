function extract_fstat(data_path, datatype_in, overwrite)

% extract_fstat: extracts the f-stat from anova data
%
% extract_fstat(data_path, datatype_in) takes the data pointed to by
% data_path (generated by path_setup) of the type datatype_in and extracts
% the f-stats from all datasets found. 
%
% extract_fstat(data_path, datatype_in, overwrite) calculates it for all
% files if overwrite is true. By default, overwrite is false, meaning that
% if the expected output file already exists, it will not be re-calculated.

% set defaults
if nargin < 3
    overwrite = false;
end

data_dir = fullfile(data_path, 'included_datasets');

% find the name of each cat
cat_names = dirsinside(data_dir);

for k = 1:numel(cat_names)
    
    % find data
    dir_sec = fullfile(data_dir, cat_names{k}, datatype_in);
    loadname = dir(fullfile(dir_sec, '*adatout.mat'));
    
    % check if the file already exists
    existname = dir(fullfile([dir_sec(1:end-1) '_fstonly/*.mat']));
    if ~overwrite && ~isempty(existname)
        fprintf('Files found for cat %s. Skipped.\n', cat_names{k})
        continue
    end
        
    for c = 1:numel(loadname)
        
        load(fullfile(dir_sec, loadname(c).name))
        
        % get data from table
        [nChan, nFreq] = size(tables);
        fstats = zeros(nChan, nFreq, 3);
        
        for ch = 1:nChan
            for f = 1:nFreq
                if isempty(tables{ch, f})
                    % table will be empty is invalid, so need to handle
                    % this error.
                    fstats(ch, f, :) = [NaN NaN NaN];
                else
                    fstats(ch, f, :) = cell2mat(tables{ch, f}(2:4, 5));
                end
            end
        end
        
        
        metavars.custom.filename = [metavars.custom.filename(1:end-4) '_fstonly'];
        
        save(metavars.custom.filename, 'metavars', 'fstats')
        
    end
    
    afpc_filemover(cat_names{k}, dir_sec, 'fstonly')

end



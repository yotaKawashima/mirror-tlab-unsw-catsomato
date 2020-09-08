function extract_fstat_y(dir_sec)

% extract_fstat: extracts the f-stat from anova data
%
% extract_fstat(data_path) takes the data pointed to by
% dir_sec and extracts the f-stats from all datasets found. 

    
% find data
loadname = dir(fullfile(dir_sec, '*adatout.mat'));

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



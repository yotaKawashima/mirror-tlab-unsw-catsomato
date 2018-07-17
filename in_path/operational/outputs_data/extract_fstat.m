function extract_fstat(data_path, datatype_in)

% directories and data types
%data_path = '/media/rannee/UNSW_Cat_Somatos/data/';
%anovadata_dir = 'collated_data/anoved_rsampsl_biprref_cmtspwr/';
%data_suffix = '_anoved_rsampsl_biprref_evkresp_cmtspwr_adatout.mat';

data_dir = fullfile(data_path, 'included_datasets');

% find the name of each cat
cat_names = dirsinside(data_dir);

for k = 1:numel(cat_names)
    
    % find data
    dir_sec = fullfile(data_dir, cat_names{k}, datatype_in, '/');
    loadname = dir(fullfile(dir_sec, '*adatout.mat'));
    
    
    for c = 1:numel(loadname)
        
        load(fullfile(dir_sec, loadname(c).name))
        
        % get data from table
        [nChan, nFreq] = size(tables);
        fstats = zeros(nChan, nFreq, 3);
        
        for ch = 1:nChan
            for f = 1:nFreq
                fstats(ch, f, :) = cell2mat(tables{ch, f}(2:4, 5));
                
            end
        end
        
        
        metavars.custom.filename = [metavars.custom.filename(1:end-4) '_fstonly'];
        
        save(metavars.custom.filename, 'metavars', 'fstats')
        
    end
    
    afpc_filemover(cat_names{k}, dir_sec, 'fstonly')

end



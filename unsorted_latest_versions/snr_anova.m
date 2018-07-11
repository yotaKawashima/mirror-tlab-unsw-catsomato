data_dir = '/media/rannee/UNSW_Cat_Somatos/data/';
out_dir = [data_dir 'collated_data/anoved_rsampsl_biprref_cmtspwr_snrsurr/'];
cd(out_dir);

% areas used in analysis
A = {'_S1', '_S2'};


% find names of cats
catnames = dirsinside([data_dir 'included_datasets/']);


% do analysis
for c = 1:numel(catnames)
    dir_sec = [data_dir 'included_datasets/' catnames{c} ...
        '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr/'];
    
    for area = 1:numel(A)
        fprintf('Calculating cat %i, area %i\n', c, area)
        filename_header = [catnames{c} A{area}];
        
        % prepare the data
        data = anp_data(dir_sec, filename_header, [0 220], false);
        
        fprintf('\tanp_data complete\n')
        
        % run the analysis
        anp_analysis(out_dir, filename_header);
        
        fprintf('\tanp_analysis complete\n')
    end
end


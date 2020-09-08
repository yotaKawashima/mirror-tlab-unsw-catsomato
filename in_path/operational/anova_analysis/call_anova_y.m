function call_anova_y(data_path, cat_name)
% Run anova2.m 
% For that, sort data first and then perform the analysis.
% This code should be able to used for log power, log SNR and VELogP.
% However, it should not be used for HGP data because the dim of HGP data
% is not the same as data of the other type due to its computation. (i.e.
% HGP is averaged across frequencies, then do not have the dim for
% frequencies.

Areas = {'_S1', '_S2'};
for a = 1:2 % per area.
    filename_header = [cat_name Areas{a}];

    % Prepare data for anova (Focus only from 0Hz to 250Hz.
    data = anp_data(data_path, filename_header, [0 250]);
    afpc_filemover(cat_name, data_path, 'adatain')

    % Performe the analysis
    data_path_a = [data_path(1:end-1), '_adatain/'];
    anp_analysis(data_path_a, filename_header);
    afpc_filemover(cat_name, data_path_a, 'adatout');
    
    % Extract only fstats
    data_path_b = [data_path_a(1:end-1), '_adatout/'];
    extract_fstat_y(data_path_b);
    afpc_filemover(cat_name, data_path_b, 'fstonly');

end

end

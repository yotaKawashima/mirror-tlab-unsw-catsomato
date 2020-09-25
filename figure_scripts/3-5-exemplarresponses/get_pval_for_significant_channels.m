% Check which bipolar channels demonstrate significance for different
% patterns. The output is saved as text file.
% Before running this script, need to run call_p_threshold

%% Data selection
% File path
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

cat_name = 'C20110808_R03'; % Section
dtypes = {'snrsurr', 'evkdpwr'}; % Data type
areas = {'S1', 'S2'}; % Area 
fois = [23, 200]; % frequencies of interest 

% Pattern
% 1=001(only interaction), 2=010(only F1), 3=011(F1 and interaction),
% 4=100(only F2), 5=101(F2 and interaction), 6=110(F1 and F2), 7=111(all)
patterns = [2, 4, 7];
patterns_type = {'only F1 main effect', 'Only F2 main effect', ...
                 'All (F1 and F2 main effect and interaction)'};

% Open Text file to write stats on
%fileID = fopen(['D:\yota\figure_scripts\4-SNRpatternS2\',...
%                'stats_significantbipolar.txt'],'w');
fileID = fopen(['./stats_significantbipolar.txt'],'w');

for i_pattern = 1:length(patterns)
    for i_dtype = 1:length(dtypes)
        for i_area = 1:length(areas)
            for i_foi = 1:length(fois)
                % Get bipolar channels showing significance for a pattern
                % for each data type, each are and each foi. 
                % Loda file
                dir_name_s = [data_path, 'collated_data/anova_props/']; 
                file_name_s = [cat_name, '_', areas{i_area},... 
                    '_F023A1F000to250_P2_',...
                    'anoved_rsampsl_biprref_evkresp_cmtspwr_',...
                    dtypes{i_dtype}, '_adatain_adatout_pthresh.mat'];
                file_s = fullfile(dir_name_s, file_name_s); 
                data_s = load(file_s);

                % Get frequency of interst
                [~, foi_ind] = find_closest(...
                                     data_s.metavars.freq{1,1}, fois(i_foi));
                
                % Get bipolar channel id 
                p_class_atfoi = data_s.p_class(:, foi_ind);
                bipolar_ch_id = find(p_class_atfoi==patterns(i_pattern));
                
                % Get p-threshold after fdr
                p_thresh = data_s.pID;                
                
                %Get statistics for the significant channel 
                % Load file
                dir_name_p = [data_path, 'included_datasets/', cat_name, ...
                    '/', 'epoched_rsampsl_biprref_evkresp_cmtspwr_',...
                    dtypes{i_dtype}, '_adatain_adatout'];
    
                file_name_p = [cat_name, '_', areas{i_area},... 
                    '_F023A1F000to250_P2_',...
                    'anoved_rsampsl_biprref_evkresp_cmtspwr_',...
                    dtypes{i_dtype},'_adatain_adatout.mat'];
                file_p = fullfile(dir_name_p, file_name_p); 
                data_p = load(file_p);

                % Get frequency of interst
                [foi_val, foi_ind] = find_closest(...
                                     data_p.metavars.freq{1,1}, fois(i_foi));

                % Get statistics (1:columns,   2:rows,       3:interaction)
                % Get statistics (1:200HZ Ampl, 2:23Hz Ampl, 3:interaction)
                pvals_atfoi = squeeze(data_p.pvals(:, foi_ind, :));
                pvals_atfoi_sig = NaN(size(pvals_atfoi));
                pvals_atfoi_sig(bipolar_ch_id,:) = ...
                                               pvals_atfoi(bipolar_ch_id,:); 
                                           
                % Write data type on the text file 
                formatSpec = ['Pattern:', patterns_type{i_pattern},...
                              ', Type:', dtypes{i_dtype}, ...
                              ', Area:', areas{i_area}, ...
                              ', Frequency of interest:%d \n'];
                fprintf(fileID, formatSpec, fois(i_foi));
                % Write down information 
                for i_bpch = 1:length(bipolar_ch_id)
                    bpch_now = bipolar_ch_id(i_bpch);
                    formatSpec = ['p-threshold:%f, bipolar ch:%3d, '...
                                  'stats:(F1 Amp,F2 Amp, interact)=',...
                                  '(%f, %f, %f) \n'];
                    A1 = [p_thresh, bpch_now, ...
                          pvals_atfoi_sig(bpch_now,2), ...
                          pvals_atfoi_sig(bpch_now,1), ...
                          pvals_atfoi_sig(bpch_now,3),];
                      
                    fprintf(fileID, formatSpec, A1);
                end % i_bpch = 1:length*biplar_ch_id
                

            end % i_foi = 1:length(fois)
        end % i_area = 1:length(areas)
    end % i_type = 1:length(types)
end % i_patterns = 1:length(patterns)

%Close text file
fclose(fileID);
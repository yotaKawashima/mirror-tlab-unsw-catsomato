% high gamma band version 3

%{
To formally analyze the high-gamma responses, we selectively removed +- 0.5
Hz around the tagged, harmonic, and intermodulatory responses in the range
of 50Hz to 150 Hz. Then, we subtracted the mean log high-gamma power in no
stimulus condition (i.e., [g,h] =[0,0]) from the log high-gamma power in
each stimulus condition, at each corresponding frequency. Finally we took
the mean across the relevant frequencies to obtain the high-gamma power.

%}

frange = [50 150];
stimf = [23 200];
k = 1;
T = 2;
% 
recordings = {'C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'};
dir_sec = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr';


nRec = numel(recordings);


%%
% for cArea = 1:2
%     
%     if cArea == 1
%         nchan = 180;
%         bipstart = 101;
%         areax = '1';
%     else
%         nchan = 112;
%         bipstart = 65;
%         areax = '2';
%     end
%     
%     
%     for cRec = 1:nRec
% 
%         [baselined, fname, chlabels, data] = hg_3d_v5_func(dir_sec, recordings{cRec}, frange, k, T, areax);
% 
%         bases{cArea}{cRec} = baselined(bipstart:end, :, :); % bipolar only
% 
%     end
%     save bases bases
% end
% 
% 
% 
% % saved in directory '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/HGP_3_tests'
% 
% % % each element in bases is a matrix: channels x trials x conditions. 
% % 
% %%
% % format for ANOVA
% for cArea = 1:2
%     
%     if cArea == 1
%         nchan = 180;
%         bipstart = 101;
%         areax = '1';
%     else
%         nchan = 112;
%         bipstart = 65;
%         areax = '2';
%     end
%     
%     for cRec = 1:nRec
% 
%         fname = dir(fullfile(dir_sec, [recordings{cRec} '*S' areax '*.mat']));
%         load(fullfile(dir_sec, fname(1).name))
% 
%         datas(cRec).fsample = data.fsample;
%         datas(cRec).custom = data.custom;
%         datas(cRec).label = data.label(bipstart:end);
% 
% 
%         y_tmp = bases{cArea}{cRec};
% 
%         d(1) = data.custom.ntrials * data.custom.subplotconfig(1);
%         d(2) = data.custom.subplotconfig(2);
%         d(3) = 1;
%         d(4) = nchan;
% 
%         % 23+trials x 200 x freq(==1) x channels
%         y = zeros(d);
%         namematrix = cell(d(1:2));
% 
%         namematrix{1, 1} = fname(1).name(24:40);
% 
% 
%         for k = 2:data.custom.conditions(2)
% 
%             col = mod(k-1, data.custom.subplotconfig(2)) +1; % column and row where the data goes
%             row = (floor((k-1)/data.custom.subplotconfig(1))*data.custom.ntrials) +1;
% 
%             dat = bases{cArea}{cRec}(:, :, k-1); % 112 x trials x 1
% 
%             y(row:row+data.custom.ntrials-1, col, 1, :) = permute(dat, [2, 4, 3, 1]);
% 
%             namematrix{row, col} = fname(k).name(24:40);
% 
%         end
% 
% 
%         datas(cRec).y = y;  
% 
% 
%         datas(cRec).conditions = namematrix;
% 
% 
%     end
%     clear data
%     save(['adatain_S' areax], 'datas')
% end

%% do anova

for cArea = 1:2
    
    clear datas
        
    if cArea == 1
        nchan = 180;
        bipstart = 101;
        areax = '1';
    else
        nchan = 112;
        bipstart = 65;
        areax = '2';
    end
    
    load(['adatain_S' areax])
    
    p_all = [];
    clear f_all
    
    for cRec = 1:nRec
        
        [pvals, tables, stats] = anp_analysis_2([], [], datas(cRec));
        mvars = rmfield(datas(cRec), 'y');
        
        save([recordings{cRec} '_S' areax '_anovahg_adatout'], 'pvals', 'tables', 'stats', 'mvars')
        
        p_all = [p_all, pvals]; % dim2 is recording number
        
        for k = 1:nchan
            for m = 1:3
                f_all(k, cRec, m) = tables{k}{m+1, 5};
            end
        end
    end
    
    save(['adatout_S' areax], 'p_all', 'f_all')
end

%% count number significant



for cArea = 1:2
    
    
        
    if cArea == 1
        nchan = 180*9;
        bipstart = 101;
        areax = '1';
    else
        nchan = 112*9;
        bipstart = 65;
        areax = '2';
    end

    clear p_all f_all
    load(['adatout_S' areax])
    
    
    [p_fdr, p_masked] = eeglab_fdr(p_all, 0.05, 'parametric');


p_class1 = double(p_masked);
p_class1(:, :, 2) = p_class1(:, :, 2)*2;
p_class1(:, :, 1) = p_class1(:, :, 1)*4;
p_class1 = sum(p_class1, 3);
    
    
    
    
    
    % standard edges. ensures no bins skipped.
    edge_st = -0.5:1:7.5;

    count_tmp = histcounts(p_class1, edge_st);

    table_values(cArea, [3, 1, 5, 2, 6, 4, 7]) = count_tmp(2:end)/nchan*100;













end




















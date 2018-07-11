for cArea = 1:2
    
    
    
    for k = 1:8
        
        filename = ['C20110808_R03_T' num2str(k, '%03i') 'T_S' num2str(cArea) '_F*A*_F*A'];
        
        data = anp_data_v2([data_dir_lvl1, 'data/HGP_data/epoched_rsampsl_biprref_cmtspwr/'], filename, [0 300], 0);
        
        
        system(['mv C20110808_R03_TStim_S' num2str(cArea) '_F000to300_P0_anov_epoched_rsampsl_biprref_cmtspwr_adatain.mat C20110808_R03_T' num2str(k, '%03i') 'T_S' num2str(cArea) '_F000to300_P0_anov_epoched_rsampsl_biprref_cmtspwr_adatain.mat']);
        
        
        
    end
    
    
    
end


%%

dir_sec = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/HGP_data/epoched_rsampsl_biprref_cmtspwr';



for k = 1%:8

    [baselined, fname, chlabels, data] = hg_3d_v4_func(dir_sec, ['C20110808_R03_T00' num2str(k) 'T'], ...
        [50 150], 1, 0.5);
    
    

end





























%%

% adatin_dir = [data_dir_lvl1, 'data/HGP_data/epoched_rsampsl_biprref_cmtspwr_adatain/'];
% fname = dir(fullfile(adatin_dir, '*.mat'));
% 
% for m = 1:numel(fname)
%     load([adatin_dir, fname(m).name])
%     
%     data.custom.filename = 
%     
% end



%%

for cArea = 1:2
    for k = 1:8
         filename = ['C20110808_R03_T' num2str(k, '%03i') 'T_S' num2str(cArea)];
         anp_analysis([data_dir_lvl1, 'data/HGP_data/epoched_rsampsl_biprref_cmtspwr_adatain/'], filename);
    end
end

%%
dir_sec = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/HGP_data/epoched_rsampsl_biprref_cmtspwr_adatout/';
filename = 'C20110808_R03_T00xT_S';
options.load = false;
options.plotchannelscan = false;
options.plotchannel = {[]; []};
options.printtype = {'-dpng'};

for k = 1:8
    filename(18) = num2str(k);
    
    [mvars, f_stats] = fstat_plot(dir_sec, filename, options, [], []);
    
    if k == 1
        ms(1).m = mvars;
        fs = f_stats;
    else
        ms(k).m = mvars;
        fs(:, k) = f_stats;
    end 
    disp(k)
end


%% extract S2 high gamma data
frange = [50 150];
stimf = [23, 200];
hars = [stimf(1):stimf(1):frange(2) stimf(2):-stimf(1):frange(1)];
hars = unique(hars);
hars = hars(hars>frange(1) & hars<frange(2));

k =0.5; 
halfrmband = (k+1)/(2*T); % remove the freq +- half the bandwidth

lband = hars - halfrmband;
hband = hars + halfrmband;

for k = 1:8
    
    % take S2 data
    f_tmp = fs{2, k};
    
    % find high gamma data
    freq_tmp = ms(k).m(2).freq{1};
    
    % find indices
    keep_i = true(1, length(freq_tmp));
    for f = 1:numel(lband)
        keep_i(freq_tmp>lband(f)&freq_tmp<hnabd(f)) = false;
    end
    
    hg(k) = mean(f_tmp(:, keep_i, :), 2);
    
    
end


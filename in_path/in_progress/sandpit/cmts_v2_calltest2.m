
filename = {'C20110511_R02';'C20110510_R05';'C20110510_R06';...
    'C20110808_R01';'C20110808_R04';'C20110808_R03';...
    'C20110808_R06';'C20110808_R09';'C20110808_Rx4'};

func_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/';
chronux_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/thirdparty_toolboxes/chronux';

addpath(genpath(func_dir))
addpath(genpath(chronux_dir))

data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';
print_format = '-dpng';

%%
input_type = '/epoched_rsampsl_biprref_evkresp/';



for m = 6%1:numel(filename)
    dir_sec = [data_dir filename{m} input_type];
    
    loadname = dir(fullfile(dir_sec, [filename{m}, '*.mat']));
    
%     params.Fs = data.fsample;
%     params.tapers = [3 5];
%     params.pad = 1;
    
    
    for k = 9%1:numel(loadname)
        load(fullfile(dir_sec, loadname(k).name))
        
        if k==1
        params.Fs = data.fsample;
        params.tapers = [3 5];
        params.pad = 1;
        end
        
        data = chronux_pwrspec_v2(data, params);
        
        
        save([data.custom.filename(1:(end-3)) 't' num2str(params.tapers(1)) '-' num2str(params.tapers(2))], 'data')
    end
    
    %afpc_filemover(filename{m}, dir_sec, 'cmtst3-5', 1)
    
end

%%
input_type = '/epoched_rsampsl_biprref_evkresp_cmtst3-5/';
secname = 'convsnr';

for m = 6%1:numel(filename)
    
    dir_sec = [data_dir filename{m} input_type];
    
    loadname = dir(fullfile(dir_sec, [filename{m}, '*.mat']));
    
    params.fmax = 300;
    params.kernalhbw = 3.2; 
    params.kernalzerolen = 1/3;
    
    
    for k = 1:numel(loadname)
        load(fullfile(dir_sec, loadname(k).name))
        
        data = conv_snr(data, params);
        
        data.custom.filename = [loadname(k).name(1:(end-4)) '_' secname];
        save(data.custom.filename, 'data')
    end
    
    afpc_filemover(filename{m}, dir_sec, secname, 1)
    
end

%%
datas = [];
input_type = '/epoched_rsampsl_biprref_evkresp_cmtst3-5_convsnr/';

for m = 1:numel(filename)
    
    dir_sec = [data_dir filename{m} input_type];
    
    loadname = dir(fullfile(dir_sec, [filename{m}, '*.mat']));   
    
    nConds = numel(loadname)/2;
    
    for k = nConds*2%*([1, 2])
        load(fullfile(dir_sec, loadname(k).name))
        disp(loadname(k).name)
        
        firstbip = prod(data.custom.spatialconfig)+1;
        
        datas = [datas; mean(data.trial(firstbip:end, :, :), 3)];
        
        
    end
    
    
    
end
filterlen = length(data.custom.convkernal);

g = mean(datas, 1, 'omitnan');
g(1:filterlen) = NaN;
g((end-filterlen):end) = NaN;
g = g(ceil(filterlen/2):(end-ceil(filterlen/2)+1));

figure(4); plot(data.freq{1}, g)
%print(gcf, print_format, 'all_cats')


%%
input_type = '/epoched_rsampsl_biprref_evkresp_cmtst3-5_convsnr/';

for m = 6
    
    dir_sec = [data_dir filename{m} input_type];
    
    loadname = dir(fullfile(dir_sec, [filename{m}, '*.mat']));   
    
    nConds = numel(loadname)/2;
    
    for k = 16
        load(fullfile(dir_sec, loadname(k).name))
        disp(loadname(k).name)
        
        firstbip = prod(data.custom.spatialconfig)+1;
        
        filterlen = length(data.custom.convkernal);
        
        data_tmp = mean(data.trial(256, :, :), 3);
        
        % remove at start/end
         data_tmp(1:filterlen) = NaN;
         data_tmp((end-filterlen):end) = NaN;
        data_tmp = data_tmp(ceil(filterlen/2):(end-ceil(filterlen/2)+1));
        
        % remove at 50, 150
        for n = [50, 150, 250]
            ind = find_closest_ind(data.freq{1}, n);
            
            data_tmp((ind-floor(filterlen/2)):(ind+floor(filterlen/2))) = NaN;
        
        end
        
        
        
        figure(6)
        plot(data.freq{1}, data_tmp)
        shg
        
    end
    
    
    
end
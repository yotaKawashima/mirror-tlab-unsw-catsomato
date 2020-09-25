% Check why some channels are invalid in S1. 
% We already know that there are not any such channels in S2 
% (plot_fstat_allfreq_yota_v4.m)
% Check what happens in time domain after bipolar re-referencing.
% Confirm that these invalid channels give 0 at the stage of 
% epoched_rsampsl_biprref

%% Set overall variables
%run(fullfile(mfilename('fullpath'), '../../path_setup.m'))
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');
datatype = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly';

%% Get channels give Nan
a = 1;
area = ['_S' num2str(a) '_'];

% find names
cat_names = dirsinside(data_dir);
k=1;
loadname = dir(fullfile(data_dir, cat_names{k}, datatype, ['*' area '*.mat']));
for k = 2:numel(cat_names)
    loadname(k) = dir(fullfile(data_dir, cat_names{k}, datatype, ['*' area '*.mat']));
end

allf = [];
sessions = [];
bpchannel_matrix = [];
% load f-stats from all conditions and store it to allf.
for k = 1:numel(loadname)
    load(fullfile(data_dir, cat_names{k}, datatype, loadname(k).name))

    % work out which channels to keep (the bipolar ones)
    if metavars.custom.nsignals == 280 || metavars.custom.nsignals == 176
        lastunipol = prod(metavars.custom.spatialconfig);
        chan = lastunipol+1:metavars.custom.nsignals;
    elseif metavars.custom.nsignals == 180 || metavars.custom.nsignals == 112
        chan = 1:metavars.custom.nsignals;
    else
        error('Only unipolar?');
    end

    allf = [allf; fstats(chan, :, :)]; % channels x freqs x 3
    sessions = [sessions ; ones(length(chan),1)*k];
    bpchannel_matrix = [bpchannel_matrix, 1:length(chan)];
end

% compute how many bipolar channels give nan value as f-stats.
% Need to check only one freq and one f-stat for each channel.
invalid_chs_concat = isnan(allf(:, 1, 1));
num_invalid_chs_concat = sum(invalid_chs_concat);
fprintf('Number of valid channels in S%i: %i out of %i\n', ... 
        a, size(allf,1) - num_invalid_chs_concat ,size(allf,1))

% Get number of invalid bipolar channels for each sessoin.
sessionid_invalidch = sessions(invalid_chs_concat);
bpchannel_invalidch = bpchannel_matrix(invalid_chs_concat);
unique_sessionids = unique(sessionid_invalidch);  

%% Load invalid channels data per session.
% Set script specific variables
file_type = 'epoched_rsampsl_biprref';
%file_type = 'epoched_rsampsl_biprref_evkresp';

for i_session = 1:length(unique_sessionids)
    session_id = unique_sessionids(i_session);
    cat_name = loadname(session_id).name(1:13);
    invalid_bpchs = bpchannel_invalidch((sessionid_invalidch == session_id));
    disp([cat_name, '-- bipolar channel id: ',...
          num2str(invalid_bpchs)]);
    
    % Find files
    loadnames = dir(fullfile(data_dir, cat_name, file_type, '*.mat'));
    nCond = numel(loadnames)/2; % number of conditions per area 
    
    % Loop through conditions in the session.
    for condid = 1:nCond
        % Load logSNR/VELogP data.
        %fprintf('Loading area %i, loading data %2i/%i\n', ...
        %    a, condid, nCond)
        load(fullfile(data_dir, cat_name, file_type, ...
            loadnames((a-1)*nCond + condid).name))
        
        % if it's the first condition, 
        % check whether data included unipolar, and
        % check if unipolar id included in srn data. If included,
        % unipolar should be removed. 
        if strcmp(data.label{1}(1:3), 'raw')
            % the first bipolar channel
            bipolarchid = 1 + prod(data.custom.spatialconfig);
            % # of bipolar channels
            nChan = data.custom.nsignals - bipolarchid + 1; 
        else
            nChan = data.custom.nsignals;
            bipolarchid = 1;
        end
        
        labels = data.label(bipolarchid:end); % label of bipolar channels
        
        % data.trial = channels x time x trials 
        data_bipolar = data.trial(bipolarchid:end, :, :);
        
        % time 
        data_time = data.time{1};
        
        % Show single channel data
        for bipolar_channel = invalid_bpchs
            % data_bipolar_now = freqs x trials
            data_bipolar_now = squeeze(data_bipolar(bipolar_channel, : ,:));
            
            % # of 0 across frequencies and trials for this channel.
            num_zero_val = sum(sum(data_bipolar_now == 0));
            
            % # of trials x # of freqs
            num_vals = size(data_bipolar_now, 1) * size(data_bipolar_now, 2);
            
            if num_zero_val ==  num_vals
                % Show nothing if all values are 0  
                %fprintf('All 0 val for %ich\n', bipolar_channel);
            elseif num_zero_val  == 0
                % No 0 value for this channel
                fprintf('No 0 val for %ich\n', bipolar_channel);
            else
                % Some 0 value for this channel
                fprintf('Some 0 val (not all) for %ich\n', bipolar_channel);
            end
            
        end % for bipolar_channel = invalid_bpchs
        %}
    end % condid=1:nCond
end % for i_session = 1:length(unique_sessionids)

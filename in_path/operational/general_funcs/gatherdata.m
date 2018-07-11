function [mixed, plotcfg] = gatherdata(cfg)

% Loads and arranges data for plotting.
%   Note: this function loads data - unlike the others in this folder.
%
%

%% configuration
nConds = length(cfg.filename);
plotcfg.nsubplots = nConds;

plotcfg.nchannels = length(cfg.channels);

% if isempty(cfg.channels)
%     cfg.channels = 1:data.custom.nsignals;
%     % if there are no channels specified, then take the last loaded value
%     % of nsignals and keep all
%     allch = true;
% end
%
% if isempty(cfg.trials)
%     cfg.trials = 1:data.custom.ntrials;
%     % if there are no trials specified, then take the last loaded value
%     % of ntrials and keep all
%     alltrials = true;
% end


%% load and condense
% preallocate for in/max
limits = zeros(nConds, 2);

for i = 1
    fprintf('\tCondition %02i\n', i)
    tic
    
    load([cfg.datadir cfg.filename{i}])
    
    mixed(i) = data;  %#ok<AGROW> no idea how many loading.
    
    allch = false;
    if isempty(cfg.channels)
        cfg.channels = 1:data.custom.nsignals;
        % if there are no channels specified, keep all
        allch = true;
    end
    
    alltrials = false;
    if isempty(cfg.trials)
        cfg.trials = 1:data.custom.ntrials;
        % if there are no trials specified, keep all
        alltrials = true;
    end
    
    % extract if needed
    if ~(allch && alltrials)
        %         mixed = rmfield(mixed, 'trial');
        
        trials = data.trial(cfg.channels, :, cfg.trials);
        
        mixed(i).trial = trials;
        
        if ~allch
            
            for j = 1:length(cfg.channels)
                chlabel{j} = data.label{cfg.channels(j)};
            end
            
            mixed(i).label = chlabel;
        end
    else
        trials = mixed(i).trial;        
    end
    
%     if ~isempty(cfg.samples)
%         trials = trials(:, cfg.samples, :);
%         mixed(i).trial = trials;
%     end
%     
    
    % remove NaNs
    trials(isinf(trials)) = NaN;
    
    limits(i, 1) = nanmin(real(trials(:)));
    limits(i, 2) = nanmax(real(trials(:)));
    taken = toc;
    fprintf('\t\tTime = %4.2f\n', taken)
end

for i = 2:nConds
    fprintf('\tCondition %02i\n', i)
    tic
    
    load([cfg.datadir cfg.filename{i}])
    
    mixed(i) = data;  %#ok<AGROW> no idea how many loading.
    
    % extract if needed
    if ~(allch && alltrials)
        %         mixed = rmfield(mixed, 'trial');
        
        trials = data.trial(cfg.channels, :, cfg.trials);
        
        mixed(i).trial = trials;
        
        if ~allch
            mixed(i).label = chlabel;
        end
        
    else
        trials = mixed(i).trial;        
    end
    
    % remove NaNs
    trials(isinf(trials)) = NaN;
    
    limits(i, 1) = min(trials(:));
    limits(i, 2) =  max(trials(:));
    
    taken = toc;
    fprintf('\t\tTime = %4.2f\n', taken)
end

limits_all(1) = min(real(limits(:, 1)));
limits_all(2) = max(real(limits(:, 2)));

plotcfg.limits = limits_all;

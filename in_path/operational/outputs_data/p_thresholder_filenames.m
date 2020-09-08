function p_thresholder_filenames(filenames, q, ifperfile)

% P_THRESHOLDER: calculates a p value threshold and returns classes
%
% p_thresholder(filenames, q) loads each file pointed to by filenames. It
%   then calculates a p-value threshold for significance using FDR. It
%   thresholds the given p-values and then saves a file with the 'classes'.
%
% p_thresholder(..., ifperfile) finds the p-value threshold individually
%   for each file when ifperfile is true. It loads all p-values and
%   determines the threshold for all files simultaneously when ifperfile is
%   false. 

if nargin < 3
    ifperfile = true;
end


pval_pool=[];

n_files = numel(filenames);

% Each file
for f = 1:n_files
    fprintf('Processing file %2i of %2i\n', f, n_files)
    % load anova result per session (load metavars)
    load(fullfile(filenames{f}))    
    % Extract only bipolar data if metavars include unipolar data.
    % pvals dim = channels x frequencies x 3 (2main effect + 1interaction)
    
    % work out which channels to keep (the bipolar ones)
    if metavars.custom.nsignals == 280 || metavars.custom.nsignals == 176
        lastunipol = prod(metavars.custom.spatialconfig);
        chan = lastunipol+1:metavars.custom.nsignals;
    elseif metavars.custom.nsignals == 180 || metavars.custom.nsignals == 112
        chan = 1:metavars.custom.nsignals;
    else
        error('Only unipolar?');
    end
       
    pvals_bipolar = pvals(chan,:,:);
    metavars.chanlabel = metavars.chanlabel(chan);
    metavars.custom.nsignals = length(chan);
    
    % FDR for each file (per sessions)
    if ifperfile
        [pID, ~] = eeglab_fdr(pvals_bipolar, q, 'parametric');
    
        % [pID,pN] = FDR(pvals,q);
        % % if pID is empty, there is no significant result
        % if isempty(pN)
        %     pN=0;
        % end
        if isempty(pID)
            pID=0;
        end
        
        pt_classifychannels(pvals, pID, q, metavars)
    % Store data for FDR all together (regardless sessions)
    else
        % Concatenate all sessions together.
        % difference frequencies, different type of f-stats, and bipolar
        % channels are all pooled together.
        all_pvals{f} = pvals_bipolar;
        pval_pool = [pval_pool;pvals_bipolar(:)];
        metas{f} = metavars;
    end
    

end

fprintf('Processing complete\n')
% FDR all together (regardless sessions)
if ~ifperfile
    % pID is a p-threshold. (scaler not vector)
    [pID, ~] = eeglab_fdr(pval_pool, q, 'parametric');
    if isempty(pID)
        pID=0;
    end
    
    for f = 1:n_files
        pt_classifychannels(all_pvals{f}, pID, q, metas{f})
    end
    
end

function pt_classifychannels(pvals, pID, q, metavars)

% first classify the channels
% p_thresh dim = channels x frequencies x 3 (2main effect + 1interaction)
p_thresh = zeros(size(pvals));
p_thresh(pvals<pID) = 1; 

% convert binary to decimal - int/(:, :, 3) is LSB
% 1=001(only interaction), 2=010(only F1), 3=011(F1 and interaction),
% 4=100(only F2), 5=101(F2 and interaction), 6=110(F1 and F2), 7=111(all)
p_class = p_thresh;
p_class(:, :, 2) = p_class(:, :, 2)*2;
p_class(:, :, 1) = p_class(:, :, 1)*(2^2);
p_class = sum(p_class, 3);

metavars.custom.filename = [metavars.custom.filename(1:end-4) '_pthresh'];

save(metavars.custom.filename, 'p_class', 'q', 'pID', 'metavars')
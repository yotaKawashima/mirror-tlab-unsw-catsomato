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


for f = 1:n_files
    fprintf('Processing file %2i of %2i\n', f, n_files)
    
    load(fullfile(filenames{f}))

    
    
    if ifperfile
        [pID, ~] = eeglab_fdr(pvals, q, 'parametric');
    
        % [pID,pN] = FDR(pvals,q);
        % % if pID is empty, there is no significant result
        % if isempty(pN)
        %     pN=0;
        % end
        if isempty(pID)
            pID=0;
        end
        
        pt_classifychannels(pvals, pID, q, metavars)
    else
        all_pvals{f} = pvals;
        pval_pool = [pval_pool;pvals(:)];
        metas{f} = metavars;
    end
    

end

fprintf('Processing complete\n')

if ~ifperfile
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
p_thresh = zeros(size(pvals));
p_thresh(pvals<pID) = 1;

% convert binary to decimal - int/(:, :, 3) is LSB
p_class = p_thresh;
p_class(:, :, 2) = p_class(:, :, 2)*2;
p_class(:, :, 1) = p_class(:, :, 1)*4;
p_class = sum(p_class, 3);

metavars.custom.filename = [metavars.custom.filename(1:end-4) '_pthresh'];

save(metavars.custom.filename, 'p_class', 'q', 'pID', 'metavars')
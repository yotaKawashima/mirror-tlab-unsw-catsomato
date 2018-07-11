function p_thresholder(data_dir, filename_out, q, ifperfile)

%loadme = '/media/rannee/UNSW_Cat_Somatos/data/included_datasets/C20110510_R05_S1_F023A1F000to250_P2_anoved_rsampsl_biprref_evkresp_cmtspwr_adatout.mat';

% data_dir = '/media/rannee/UNSW_Cat_Somatos/data/included_datasets/';
% a = 1; cat_name = 'C20110510_R05';
% filename_out = [cat_name '*S' num2str(a)];
% q=0.05;




if nargin < 4
    ifperfile = true;
end

loadname = dir(fullfile(data_dir, [filename_out '*adatout.mat']));


pval_pool=[];


for f = 1:numel(loadname)
    load(fullfile(data_dir, loadname(f).name))

    
    
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

if ~ifperfile
    [pID, ~] = eeglab_fdr(pval_pool, q, 'parametric');
    if isempty(pID)
        pID=0;
    end
    
    for f = 1:numel(loadname)
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
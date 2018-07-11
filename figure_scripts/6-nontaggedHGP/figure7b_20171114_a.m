%% load

for a = 1:2
A = num2str(a);
load(['C20110808_R03_S' A '_FxxxAxxx_FxxxAxxx_epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp.mat'])

%% preproc
% work out which channels to keep (the bipolar ones)
if strcmp(data.label{1}(1:3), 'raw')
    c1_ind = 1 + prod(data.custom.spatialconfig);
    nChan = data.custom.nsignals - c1_ind + 1;
else
    nChan = data.custom.nsignals;
    c1_ind = 1;
end

% frequency is already selected

labels{a} = data.label(c1_ind:end);
spatialconfigs{a} = data.custom.spatialconfig;
condfigs{a} = data.custom.subplotconfig;

%% loop through conditions
for t = 1:numel(data.trial)
    %% select the right channels
    hgp_power = mean(data.trial{t}(c1_ind:end, :, :), 2);
    
    %% perform t-test
    for k = 1:nChan
        [allhs{a}(k,1,t) , allps{a}(k,1,t)]=ttest(hgp_power(k, 1, :), 0, 'Tail', 'right');
        
        
    end
    
    %% find means
    mean_hgp{a}(:, :, t) = mean(hgp_power, 3);
end



%% threshold
q = 0.05;


mean_hgp{a}(allhs{a}==0)=NaN;


% do FDR
[pID, ~] = eeglab_fdr(allps{a}(:, 1, :), q, 'parametric');
tmp = mean_hgp{a}(:, 1, :);
tmp(allps{a}(:, 1, :)>=pID)=NaN;

mean_hgp{a}(:, 1, :) = tmp;



%% plot
fnum=1; % first figure number
cmap = hot(64);
cmap = cmap(1:60, :);

nAreas = numel(mean_hgp);
[nChans, ~, nConds] = size(mean_hgp{1});


clim = find_lims(mean_hgp);
clim(2) = clim(2);

    % relabel channels
    relabels = draw_biprref_chlabfunction(labels{a});
    condfig = condfigs{a};
    

        fig = (a-1)*nAreas + fnum; % has bugs if too high but works for now
        % make a figure
        figure(fig)
        clf
        set(gcf, 'Name', ['mean HGP S' num2str(a)]) 
        
        for co = 1:nConds
            % make subplot
            subtightplot(condfig(1), condfig(2), co+1)
            % goes across rows
            
            draw_biprref(mean_hgp{a}(:,1,co), relabels, spatialconfigs{a}, clim)
            colormap(cmap)
        end
        
        % call function to clean up and label the plot
        subtightplotcleaner(fig, condfig, 'cleanticks', false, ...
            'catnames', [{'F023A000_F200A000'}, data.datalabels], ...
            'topinds', 10:17, 'sideinds', 1:8, 'box', false)
        
        dt = datestr(now, 'yyyymmdd');
        
        % save figure
         print(gcf, '-depsc', ['Figure7_S' num2str(a) '_' dt])


end

figure((a-1)*nAreas +2)
colorbar
colormap(cmap)
caxis(clim)
print(gcf, '-depsc', ['Figure7_colorbar_' dt])

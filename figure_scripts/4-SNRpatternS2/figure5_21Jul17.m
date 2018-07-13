

cat = 'C20110808_R03';
FOI = [23, 200];

data_dir = '/media/rannee/UNSW_Cat_Somatos/data/included_datasets/';
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';

% find list of names
loadnames = dir(fullfile(data_dir, cat, data_type, '*.mat'));

% find some numbers
nCond = numel(loadnames)/2;
nFOI = numel(FOI);

allps = cell(1, 2); % one for each area
allhs = cell(1, 2);
meansnrs = cell(1, 2);

for a = 1:2
    for c = 1:nCond
        %disp((a-1)*nCond + c)
        
        
        % load
        load(fullfile(data_dir, cat, data_type, loadnames((a-1)*nCond + c).name))
        
        % if it's the first, preallocate some things
        if c==1
            % check if unipolar
            if strcmp(data.label{1}(1:3), 'raw')
                c1_ind = 1 + prod(data.custom.spatialconfig);
                nChan = data.custom.nsignals - c1_ind + 1;
            else
                nChan = data.custom.nsignals;
                c1_ind = 1;
            end
            allps{a} = zeros(nChan, nFOI, nCond);
            allhs{a} = zeros(nChan, nFOI, nCond);
            meansnrs{a} = zeros(nChan, nFOI, nCond);
            
            iFOI = zeros(1, nFOI);
            vFOI = zeros(1, nFOI);
            for k = 1:nFOI
                [vFOI(k), iFOI(k)] = find_closest(data.freq{1}, FOI(k));
            end
            
        else
            % check that the value of frequency is right
            if ~isequal(data.freq{1}(iFOI), vFOI')
                warning('Frequency indicies inconsistent. Recalculating. (a=%i, c=%i)', a, c)
                
                for k = 1:nFOI
                    [vFOI(k), iFOI(k)] = find_closest(data.freq{1}, FOI(k));
                end
            end
        end
        
        % extract snrs
        snrs = data.trial(c1_ind:end, iFOI, :);
                
        % do the t-test
        for k = 1:nChan
            for f = 1:nFOI
                [allhs{a}(k,f,c) , allps{a}(k,f,c)]=ttest(snrs(k, f, :), 0, 'Tail', 'right');
                
            end
        end
        
        % find mean
        meansnrs{a}(:, :, c) = mean(snrs, 3);
    end
    
    labels{a} = data.label(c1_ind:end);
    spatialconfigs{a} = data.custom.spatialconfig;
    condfigs{a} = data.custom.subplotconfig;
end

save f5dat allhs allps meansnrs labels spatialconfigs condfigs


%% threshold
q = 0.05;

for a = 1:2
    meansnrs{a}(allhs{a}==0)=NaN;
    
    for f = 1:nFOI
        % do FDR
        [pID, ~] = eeglab_fdr(allps{a}(:, f, :), q, 'parametric');
        tmp = meansnrs{a}(:, f, :);
        tmp(allps{a}(:, f, :)>=pID)=NaN;
        
        meansnrs{a}(:, f, :) = tmp;
    end
    
end

%% visualise

fnum=1; % first figure number
cmap = 'hot';

nAreas = numel(meansnrs);
[nChans, nFOI, nConds] = size(meansnrs{1});

clim = find_lims(meansnrs);

for a = 1:nAreas
    % relabel channels
    relabels = draw_biprref_chlabfunction(labels{a});
    condfig = condfigs{a};
    
    for f = 1:nFOI
        fig = (a-1)*nAreas + (f-1) + fnum; % has bugs if too high but works for now
        % make a figure
        figure(fig)
        clf
        set(gcf, 'Name', ['f' num2str(FOI(f)) ' S' num2str(a)]) 
        
        for co = 1:nConds
            % make subplot
            subtightplot(condfig(1), condfig(2), co)
            % goes across rows
            
            draw_biprref(meansnrs{a}(:,f,co), relabels, spatialconfigs{a}, clim)
            colormap(cmap)
        end
        
        % call function to clean up and label the plot
        subtightplotcleaner(fig, condfig, 'cleanticks', false, ...
            'catnames', {loadnames((a-1)*nConds+1:a*nConds).name}, ...
            'topinds', 27:34, 'sideinds', 18:25, 'box', false)
        
        % save figure
        print(gcf, '-depsc', ['Figure5_f' num2str(FOI(f)) '_S' num2str(a)])
    end
end



figure((a-1)*nAreas + (f-1) + fnum+1)
colorbar
colormap(cmap)
caxis(clim)
print(gcf, '-depsc', 'Figure5_colorbar')
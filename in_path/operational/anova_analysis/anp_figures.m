function anp_figures(data_dir, filename_header, options)
% anp_figure: draws figures from data in the anova2 pipeline
%   anp_figures(data_dir, filename_header, options) draws the figures
%   specified by the struct options for the data in data_dir, specified by
%   filename_header. 
%   
%   The type of plot produced depends on the struct options. See code for
%   more details. For multiple plots, simply use an array of structures.
%
%   Currently if no classification limit for the p-values is provided, a
%   FDR with Q = 1% is used instead.

% Open file anp_figureshelper.m for templates of options struct. Otherwise,
% it generally has the following fields: 
%   (1x1 double) clim, classlim
%   (1x1 logical) singfoi, singchan, classify
%   (row vector) foi, chan

% Written by Rannee Li, Jul 2015

% maybe able to replace eval statements with .() notation? 

%% load data
fname = dir(fullfile(data_dir, [filename_header '*_P*_anov*_adatout.mat']));
load([data_dir fname.name])

%% process metavars
[frow, fcol] = a2c_depfreq(metavars);
% allocate the dependence names
dep_names = {[fcol ' Hz dependence'], [frow ' Hz dependence'], 'interaction'};

%% find the correct subfunction and do some preprocessing
for rep = 1:numel(options)
    sel = bi2de([options(rep).singfoi options(rep).singchan options(rep).classify], 'left-msb');
    
    switch sel
        case 0 % all frequencies, all channels, p-values plot
            switch options(rep).clim
                case 'default'
                    options(rep).clim = [];
                case 'fdr'
                    q = options(rep).fdr{1};
                    [pID,pN] = FDR(pvals,q); %#ok<ASGLU,NASGU> may be used in eval below
            
                    eval(['options(rep).clim = ' options(rep).fdr{2} ';']);
            end
            
            anpf_FA_CHA_noclass(pvals, metavars, options(rep), dep_names, fname.name)
        case 1 % all frequencies, all channels, classify by p-value
            % find out if FDR is needed
            if numel(options(rep).fdr) > 1 % then fdr is true
                q = options(rep).fdr{1};
                [pID,pN] = FDR(pvals,q); %#ok<ASGLU,NASGU> may be used in eval below
                % if pID is empty, there is no significant result 
                % Addded by Phoebe You 20Jan2016
                if isempty(pID), pID=0; end
                if isempty(pN), pN=0; end             
                eval(['options(rep).classlim = ' options(rep).fdr{2} ';']);
            end
            
            anpf_FA_CHA_yeclass(pvals, metavars, options(rep), frow, fcol, fname.name)
        case 2 % all frequencies, one channel, p-values plot
            % check the type of clim values 
            if isnumeric(options(rep).clim)
                % assume it is 'default'
                options(rep).clim = [];
            else
                if isfield(options(rep), 'fdr') && numel(options(rep).fdr)>1
                    % use fdr to determine one of the horizontal lines
                    q = options(rep).fdr{1};
                    [pID,pN] = FDR(pvals,q); %#ok<ASGLU,NASGU> may be used in eval below
            
                    eval(['p_calc = ' options(rep).fdr{2} ';']);
                    
                    options(rep).clim = [options(rep).clim p_calc];
                end
            end
            
            anfp_FA_CH1_noclass(pvals, metavars, options(rep), dep_names, fname.name)
        case 4 % one frequency, all channels, p-values plot
            switch options(rep).clim
                case 'default'
                    options(rep).clim = [];
                case 'fdr'
                    q = options(rep).fdr{1};
                    [pID,pN] = FDR(pvals,q); %#ok<ASGLU,NASGU> may be used in eval below
            
                    eval(['options(rep).clim = ' options(rep).fdr{2} ';']);
            end
            
            anfp_F1_CHA_noclass(pvals, metavars, options(rep), dep_names, fname.name)
        case 5 % one frequency, all channels, classify by p-value
            % find out if FDR is needed
            if numel(options(rep).fdr) > 1 % then fdr is true
                q = options(rep).fdr{1};
                [pID,pN] = FDR(pvals,q); %#ok<ASGLU,NASGU> may be used in eval below
                
                eval(['options(rep).classlim = ' options(rep).fdr{2} ';']);
            end
            
            anfp_F1_CHA_yeclass(pvals, metavars, options(rep), frow, fcol, fname.name)
        otherwise
            fprintf('You selected sel = %i.\n', sel)
            error('Desired type of plot is meaningless.')
    end
end

%% Subfunction: sel == 0
function anpf_FA_CHA_noclass(pvals, metavars, options, dep_names, imgname)
% Plots p-value (colour) for response frequency (x) against channel (y)
%   This image shows all channels and all frequencies

% if there are no limits provided, then use the default.
if ~isfield(options, 'clim') || isempty(options.clim)
    clim = [];
else
    clim = [0, options.clim];
end

for l = 1:3
figure(l); clf
imagesc(metavars.freq{1},[], pvals(:, :, l), clim)
% imagesc(metavars.freq{1}, pvals(1, :, l), clim)
colorbar
colormap default


title({'p-value for all channels across all frequencies', dep_names{l}})
ylabel('channels')
xlabel('response frequency (Hz)')

print(l, '-dpng', [imgname(1:end-4) '_CallFall_pval' dep_names{l}(1:3)])
end

%% subfunction: sel == 1
% this subfunction has been moved to a separate file for high gama analysis

% function anpf_FA_CHA_yeclass(pvals, metavars, options, frow, fcol, imgname)
% % Plots classification (colour) for response fq (x) against channel (y)
% %   This image shows all channels and all frequencies
% 
% % first classify the channels
% p_thresh = zeros(size(pvals));
% p_thresh(pvals<options.classlim) = 1;
% % convert binary to decimal - int/(:, :, 3) is LSB
% p_class = p_thresh;
% p_class(:, :, 2) = p_class(:, :, 2)*2;
% p_class(:, :, 1) = p_class(:, :, 1)*4;
% p_class = sum(p_class, 3);
% 
% figure(1); clf
% imagesc(metavars.freq{1}, [], p_class)
% 
% colorbar
% tab_base = de2bi(0:7, 'left-msb');
% colormap(tab_base)
% colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
% [fcol '+int'], [frow '+' fcol], 'all'})
% 
% title({'Classification of channel by p-value across all frequencies', ['p < ' num2str(options.classlim)]})
% ylabel('channels')
% xlabel('response frequency (Hz)')
% 
% p_str = num2str(options.classlim);
% p_str(2) = 'p';
% print(1, '-dpng', [imgname(1:end-4) '_CallFall_plt' p_str])

%% subfunction: sel == 2
% function anfp_FA_CH1_noclass(pvals, metavars, options, dep_names, imgname)
% % Plots response frequency (x) against p-value (y)
% %   This image shows one channel and all frequencies
% fq = metavars.freq{1};
% log_p = log10(pvals);
% 
% for ch = 1:numel(options.chan) % loop across all channels
%     figure(ch); clf
%     hold on
%     
%     % plot data
%     ch_ind = options.chan(ch);
%     plot(fq, log_p(ch_ind, :, 1), 'r', fq, log_p(ch_ind, :, 2), 'g', ...
%         fq, log_p(ch_ind, :, 3), 'b')
%     legend(dep_names, 'Location', 'southoutside', 'Orientation', 'horizontal')
%     
%     % draw foi lines
%     if isfield(options, 'foi') && ~isempty(options.foi)
%         screenheight = get(gca, 'YLim');
%         for f = 1:numel(options.foi)
%             hold on
%             plot(options.foi(f)*[1, 1], screenheight, 'k-.')
%             hold off
%         end
%     end
%     
%     % draw p-value lines
%     if isfield(options, 'clim') && ~isempty(options.clim)
%         screenwidth = get(gca, 'XLim');
%         for f = 1:numel(options.clim)
%             hold on
%             plot(screenwidth, log10(options.clim(f))*[1, 1], 'k-.')
%             hold off
%         end
%     end
%       
%     
%     hold off
%     
%     ylabel('log10(p)')
%     xlabel('frequency (Hz)')
%     title({'p-value for all frequencies'; ['Ch: ' metavars.chanlabel{ch_ind}]}, 'interpret', 'none')
% 
%     print(ch, '-dpng', [imgname(1:end-4) '_C' num2str(ch_ind, '%03i') 'Fall_pval' dep_names{ch}(1:3)])
% end

%% subfunction: sel == 4
function anfp_F1_CHA_noclass(pvals, metavars, options, dep_names, imgname)
% Plots p-values by channel arrangement. 
%   This image shows one channel and one frequency

ch_label = draw_biprref_chlabfunction(metavars.chanlabel);
if isfield(options, 'clim') && ~isempty(options.clim)
    clim = [0, options.clim];
else
    clim = [];
end

for k = 1:numel(options.foi)
    [~, foi_ind] = find_closest(metavars.freq{1}, options.foi(k));
    for l = 1:3
        figure(l); clf
        draw_biprref(pvals(:, foi_ind, l), ch_label, metavars.custom.spatialconfig, clim)
        title(['p-values for response dependance on ' dep_names{l} ' at ' num2str(options.foi(k)) ' Hz'])
        colorbar
        colormap default
        
        print(l, '-dpng', [imgname(1:end-4) '_CallF' num2str(options.foi(k)) '_pval' dep_names{l}(1:3)])
    end
    
    
end

%% subfunction: sel == 5
function anfp_F1_CHA_yeclass(pvals, metavars, options, frow, fcol, imgname)
% Classifies channel in channel arrangement. 
%   This image shows one channel and one frequency

ch_label = draw_biprref_chlabfunction(metavars.chanlabel);

% first classify the channels
p_thresh = zeros(size(pvals));
p_thresh(pvals<options.classlim) = 1;
% convert binary to decimal - int/(:, :, 3) is LSB
p_class = p_thresh;
p_class(:, :, 2) = p_class(:, :, 2)*2;
p_class(:, :, 1) = p_class(:, :, 1)*4;
p_class = sum(p_class, 3);

for k = 1:numel(options.foi)
    [~, foi_ind] = find_closest(metavars.freq{1}, options.foi(k));
    
    figure(k); clf
    draw_biprref(p_class(:, foi_ind), ch_label, metavars.custom.spatialconfig)
    title({'Classification of channel by ANOVA2 results', ['(p < ' num2str(options.classlim) ')']})

    colorbar
    tab_base = de2bi(0:7, 'left-msb');
    colormap(tab_base)
    colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
        [fcol '+int'], [frow '+' fcol], 'all'})
    
    p_str = num2str(options.classlim);
    p_str(2) = 'p';
    print(k, '-dpng', [imgname(1:end-4) '_CallF' num2str(options.foi(k)) '_plt' p_str])
end



%%
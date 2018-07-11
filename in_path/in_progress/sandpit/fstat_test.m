addpath(genpath('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/'))

data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';
cat = 'C20110808_R03';

ifload = false; % true; %
ifplot = true;

% plot with f-statistic
fname = dir(fullfile(data_dir, data_type, [cat '*'])); % should return 2 files: S1 and S2



if ifload
    f_stats = cell(2, 1);
    
    for k = 1:2 % loop through S1 and S2
        load(fullfile(data_dir, data_type, fname(k).name))

        % The f-statistic is stored in tables, which are massive cell arrays. This
        % should be fun.

        % preallocate 
        [d1, d2] = size(tables);
        f_stats{k} = zeros([d1, d2, 3]); % always 3 values: 23, 200, int

        % LOOOOOOP because cell arrays
        for m = 1:d1
            for n = 1:d2
                for p = 1:3
                    f_stats{k}(m, n, p) = tables{m, n}{p+1, 5};
                end
            end
        end

        mvars(k) = metavars;
        
    end
    
    savename = fname(1).name;
    savename(22) = 'x';
    savename = [savename(1:end-4) '_2fstats'];
    save([fname(1).name], 'f_stats', 'mvars')
    
end

if ifplot
    order = {'200', '023', 'int'};
    
    for k = 1:2
        printname1 = fname(k).name(1:22);
        printname = [printname1 '_fstatplot'];        
        
        clim = find_lims(f_stats{k});

        for p = 1:3
            figure((k-1)*3+p); clf
            imagesc(mvars(k).freq{1}, [], f_stats{k}(:, :, p), clim)
            colorbar
            colormap jet

            xlabel('frequency (Hz)')
            ylabel('channels')
            
            print(gcf, '-dpng', [printname '_' order{p}])

        end

    end
end



figure(8);plot(mvars(2).freq{1}, permute(f_stats{2}(101, :, :), [2, 3, 1]))
xlabel('frequency (Hz)')
ylabel('f-statistic')
title({'F-statistic from ANOVA analysis'; [fname(1).name(1:13) ', ' fname(1).name(21:22) ', channel ' mvars(1).chanlabel{101}]}, 'Interpreter', 'none')
legend(order, 'Location', 'SouthOutside', 'Orientation', 'horizontal')
print(gcf, '-dpng', [printname '_CH165'])
function [mvars, f_stats] = fstat_plot(data_dir, filename, options, mvars, f_stats)

% Extracts and plots f-statistic from computed ANOVA data.
% 
%   fstat_plot(data_dir, filename, options) loads the data from the
%   corresponding to filename (string) in the directory specified by the
%   string data_dir. data_dir is the full data directory, i.e. the file is
%   directly inside data_dir. Options controls the types of plots produced
%   and the print format. See below for more details on struct options.
%
%   fstat_plot(..., mvars, fstats) loads data if either or both of mvars
%   and fstats are empty, otherwise it uses the data passed in. data_dir
%   and filename are needed regardless of whether the files are to be
%   loaded, as they are used to name the output files.
% 
%   [mvars, f_stats] = fstat_plot(...) outputs the data.
%
% Fields in options
%   load: 
%       1x1 logical. If true, the data is retrieved from data_dir and
%       loaded. If false, this value is overwritten if no data is passed to
%       the function.
%   plotchannelscan:
%       1x1 logical. If true, plots results for all channels using imagesc. 
%   plotchannel:
%       2 element cell. Element 1 contains S1 channels; Element 2 contains
%       S2 channels. Each cell element contains a vector containing
%       integers (doubles) that specify which single channel plots to
%       produce. Otherwise, leave as []
%   printtype: 
%       '-dpng', '-epsc', or any combination of any flags used in the print
%       function.

% Written by Rannee Li, Jul 2016

% find the filenames. Even if we don't load them, we still need them to
% write files to output using the names of the input
fname = dir(fullfile(data_dir, [filename '*']));

if nargin<4 || numel(mvars)<1 || numel(f_stats)<1 || options.load
    % if we need to load data, basically
    f_stats = cell(2, 1);
    
    for k = 1:2 % loop through S1 and S2
        load(fullfile(data_dir, fname(k).name))

        % The f-statistic is stored in tables, which are massive cell
        % arrays. This should be fun.

        % preallocate 
        [d1, d2] = size(tables); %#ok<USENS> variable is loaded
        f_stats{k} = zeros([d1, d2, 3]); % always 3 values: 23, 200, int

        % LOOOOOOP because cell arrays
        for m = 1:d1
            for n = 1:d2
                for p = 1:3
                    f_stats{k}(m, n, p) = tables{m, n}{p+1, 5};
                end
            end
        end
        
        if k == 1
            % push metavars to output because we need them for graphing
            mvars = metavars;
        else
            mvars(2) = metavars;
        end
        
    end
    
    % save the files to output so we don't lose our precious calculations
    savename = fname(1).name;
    savename(22) = 'x';
    savename = [savename(1:end-4) '_2fstats'];
    save(savename, 'f_stats', 'mvars')
end

nameorder = {'200', '023', 'int'}; % variable used for naming the printouts

% plot channelscan
if options.plotchannelscan
    for k = 1:2
        % name for the printout
        printname1 = fname(k).name(1:22);
        printname = [printname1 '_fstatplot'];

        clim = find_lims(f_stats{k});

        for p = 1:3
            figure((k-1)*3+p); clf
            imagesc(mvars(k).freq{1}, [], f_stats{k}(:, :, p), clim)
            colorbar
            colormap hot

            xlabel('frequency (Hz)')
            ylabel('channels')

            for m = 1:numel(options.printtype)
                print(gcf, options.printtype{m}, [printname '_' nameorder{p}])
            end

        end

    end
end

% plot channels
for k = 1:2
    for n = 1:numel(options.plotchannel{k})
        figure()

        plot(mvars(n).freq{1}, permute(f_stats{k}(options.plotchannel{k}(n), :, :), [2, 3, 1]))

        xlabel('frequency (Hz)')
        ylabel('f-statistic')
        title({'F-statistic from ANOVA analysis'; [fname(k).name(1:13) ', ' fname(k).name(21:22) ', channel ' mvars(k).chanlabel{options.plotchannel{k}(n)}]}, 'Interpreter', 'none')
        legend({'200 Hz', '23 Hz', 'interaction'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')

        for m = 1:numel(options.printtype)
            print(gcf, options.printtype{m}, [printname '_CH' num2str(options.plotchannel{k}(n), '%03i')])
        end
    end
end

if nargout < 2
    clear mvars f_stats
end










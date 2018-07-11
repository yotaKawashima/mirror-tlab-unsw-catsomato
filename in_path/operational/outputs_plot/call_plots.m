function call_plots(mixed, plotcfg, plotmethod, filename)
% plots using subplotconfig in data.custom

set(0,'DefaultTextInterpreter','none');
% For some reason, the yaxis direction is reversed.

% fprintf('conditioning\n')
if strcmp(plotmethod, 'plot')
    ifplot = 1;
    colors = {'b', 'k', 'r', 'y', 'm', 'c', 'g'};
elseif strcmp(plotmethod, 'imagesc')
    ifplot = 2;
elseif strcmp(plotmethod, 'biprref')
    ifplot = 3;
elseif strcmp(plotmethod, 'array')
    ifplot = 4;
else
    error('input plotmethod invalid')
end

% fprintf('xaxis\n')
try
    xaxis = mixed(1).time{1};
catch %#ok
    xaxis = mixed(1).freq{1};
end


% fprintf('figure\n')
for i = 1:mixed(1).custom.conditions(2)
%     fprintf('\tdrawing %02i\n', i)
    subtightplot(mixed(1).custom.subplotconfig(2), mixed(1).custom.subplotconfig(1), i)
    
    hold on
    
    if ifplot == 1
        for j = 1:plotcfg.nchannels
            figure(j)
            subtightplot(mixed(1).custom.subplotconfig(2), mixed(1).custom.subplotconfig(1), i)
            plot(xaxis, squeeze(mean(mixed(i).trial(j, :, :), 3)), colors{mod(j, 7) +1})
            set(gca, 'ylim', [plotcfg.limits(1) plotcfg.limits(2)])
            
        end
    elseif ifplot == 2
        imagesc(xaxis, [], mean(real(mixed(i).trial), 3), [plotcfg.limits(1) plotcfg.limits(2)])
    elseif ifplot == 3
        draw_biprref_v2(mixed(i), plotcfg.biprref_index, plotcfg)
    elseif ifplot == 4
         draw_array(mixed(i), plotcfg)
    end
    hold off
    
    axis tight
    
    
    if i == (mixed(1).custom.conditions(2) - mixed(1).custom.subplotconfig(1) + 1)
        set(gca, 'FontSize', 6)
        
%        ylabel('ylabel-TBA')
%        xlabel('xlabel-TBA')
        
%         xlabh = get(gca,'XLabel'); set(xlabh,'Position',get(xlabh,'Position') + [0 5 0]) % move close for print
        ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [0.2 0 0]) % move close for print
    else
        set(gca, 'XTickLabel', [], 'YTickLabel', [])
    end
    
    
    if i <= mixed(1).custom.subplotconfig(1)
        title(mixed(i).custom.filename(33:40))
    end
    
    
    if mod(i, mixed(1).custom.subplotconfig(1))==0
        set(gca, 'yaxislocation', 'right')
        ylabel(mixed(i).custom.filename(24:31))
    end
    
end

if nargin > 3 % if filename is an input argument
    print(gcf, '-dpng', filename)
end
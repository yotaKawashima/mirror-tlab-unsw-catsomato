% collapses the fstat correlation into "bars"

if verLessThan('matlab', '8.4')
    error('Code does not work as errorbar cannot be plotted horizontally.')
elseif verLessThan('matlab', '9.2')
    warning('Code may not work if errorbar cannot be plotted horizontally/')
end

figure(1); clf;

for a = 1:2
    area = num2str(a);
    imgformat = '-dpng';
    
    
    load(['S' area '_corrdata.mat']);
    
    [R, P, RL, RU] = corrcoef(log(corrdata));
    
    % preallocate
    x = zeros(3, 10);
    y = repmat(1:10, 3, 1)+(a-1)*0.2;
    xl = zeros(3, 10);
    xu = zeros(3, 10);
    
    % row 1: 23 to 200.
    t = 1;
    r1 = 1;
    c1 = 11;
    x(:, t)  = diag( R(r1:20:end, c1:20:end));
    xl(:, t) = diag(RL(r1:20:end, c1:20:end));
    xu(:, t) = diag(RU(r1:20:end, c1:20:end));
    
    % row 2: 23 to harmonics
    t = 2;
    r1 = 1;
    c1 = 2:10;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    % row 3: 23 to intermodulation
    t = 3;
    r1 = 1;
    c1 = 12:19;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    % row 4: 23 to HGP
    t = 4;
    r1 = 1;
    c1 = 20;
    x(:, t)  = diag( R(r1:20:end, c1:20:end));
    xl(:, t) = diag(RL(r1:20:end, c1:20:end));
    xu(:, t) = diag(RU(r1:20:end, c1:20:end));
    
    % row 5: 200 to harmonics
    t = 5;
    r1 = 11;
    c1 = 2:10;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    % row 6: 200 to intermodulation
    t = 6;
    r1 = 11;
    c1 = 12:19;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    % row 7: 200 to HGP
    t = 7;
    r1 = 11;
    c1 = 20;
    x(:, t)  = diag( R(r1:20:end, c1:20:end));
    xl(:, t) = diag(RL(r1:20:end, c1:20:end));
    xu(:, t) = diag(RU(r1:20:end, c1:20:end));
    
    % row 8: harmonics to intermodulation
    t = 8;
    r1 = 12:19;
    c1 = 2:10;
    for k = 1:3
        x(k, t)  = mean(mean( R(r1+(k-1)*20, c1+(k-1)*20)));
        xl(k, t) = mean(mean(RL(r1+(k-1)*20, c1+(k-1)*20)));
        xu(k, t) = mean(mean(RU(r1+(k-1)*20, c1+(k-1)*20)));
    end
    
    % row 9: harmonics to HGP
    t = 9;
    r1 = 2:10;
    c1 = 20;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    % row 10: intermodulation to HGP
    t = 10;
    r1 = 12:19;
    c1 = 20;
    for k = 1:3
        x(k, t)  = mean( R(r1+(k-1)*20, c1+(k-1)*20));
        xl(k, t) = mean(RL(r1+(k-1)*20, c1+(k-1)*20));
        xu(k, t) = mean(RU(r1+(k-1)*20, c1+(k-1)*20));
    end
    
    
    xl = x - xl;
    xu = xu - x;
    
    
    % plot
    
    mark = 'os';
    col = 'rgb';
    hold on
    for k = 1:3
        h = errorbar(x(k, :), y(k, :), xl(k, :), xu(k, :), 'horizontal', mark(a));
        set(h, 'LineWidth', 2, 'Color', col(k))
    end
    legend({'Main Effect 23', 'Main Effect 200', 'Interaction'}, ...
        'Location', 'SouthOutside', 'Orientation', 'Horizontal')
    hold off
    
%     ax(a) = gca;
%     
%     xlim_tmp(a, :) = get(gca, 'XLim');
    
    set(gca, 'YDir', 'reverse');
end

title('S1 = circle, S2 = square')

ytl = {'1) f1-f2', '2) f1-harm' , '3) f1-IM', '4) f1-HGP', ...
    '5) f2-harm', '6) f2-IM', '7) f2-HGP', '8) harm-IM', '9) harm-HGP', ...
    '10) IM-HGP', ''};

set(gca, 'YTick', 1.1:11.1, 'YTickLabel', ytl, 'YLim', [0.5, 10.7])

print(gcf, imgformat, ['S1and2_correlationbars']);
% xlim(1) = min(xlim_tmp(:, 1));
% xlim(2) = max(xlim_tmp(:, 2));
% 
% 
% for a = 1:2
%     figure(a)
%     set(ax(a), 'XLim', xlim);
%     title(['S' num2str(a)])
%     print(gcf, imgformat, ['S' num2str(a) '_correlationbars']);
%     
% end
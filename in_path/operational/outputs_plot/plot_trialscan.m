function plot_trialscan(data_dir, erpname, pwrname, channel)

% figure 1 for erp, figure 2 for power
CondNames = {'F023A000_F200A016', 'F023A159_F200A000', 'F023A159_F200A016'};

plotdat(3).erpdat = []; % preallocate

% load and extract data
for i = 1:3
    condname = CondNames{i};
    
    % load ERP
    data = trialscan_load(data_dir, erpname, condname);
    [xax, twod] = trialscan_data(data, channel);
    plotdat(i).erpdat = twod;
    plotdat(i).erpxax = xax;
    
    % load power
    data = trialscan_load(data_dir, pwrname, condname);
    [xax, twod] = trialscan_data(data, channel);
    plotdat(i).pwrdat = twod;
    plotdat(i).pwrxax = xax;
end

clear data

% fimd limits for plotting
erplims = find_lims(plotdat, 'erpdat');
pwrlims = find_lims(plotdat, 'pwrdat');

figure(1); clf
figure(2); clf

% plot
for i = 1:3
    figure(1)
    subplot(2, 2, i+1)
    imagesc(plotdat(i).erpxax, [], plotdat(i).erpdat, erplims)
    title(CondNames{i}, 'Interpret', 'none')
    ylabel('trial')
    xlabel('time (sec)')
    
    figure(2)
    subplot(2, 2, i+1)
    imagesc(plotdat(i).pwrxax, [], plotdat(i).pwrdat, pwrlims)
    title(CondNames{i}, 'Interpret', 'none')
    ylabel('trial')
    xlabel('frequency (Hz)')
end

figure(1)
suptitle(['ERPs by trial for channel ' num2str(channel, '%03i') ' in ' erpname(21:22)])
figure(2)
suptitle(['Power spectrum by trial for channel ' num2str(channel, '%03i') ' in ' erpname(21:22)])

% do colorbars
figure(3); clf; 
colorbar
caxis(erplims)
title(['Colorbar for: ERPs by trial for channel ' num2str(channel, '%03i') ' in ' erpname(21:22)])

figure(4); clf; 
colorbar
caxis(pwrlims)
title(['Colorbar for: Power spectrum by trial for channel ' num2str(channel, '%03i') ' in ' erpname(21:22)])




% % start plotting!
% % 200 only
% condname = 'F023A000_F200A016';
% 
% figure(1)
% subplot(2, 2, 2)
% lims = trialscan_helper(data_dir, erpname, condname, channel);
% erplim = lims;
% 
% figure(2)
% subplot(2, 2, 2)
% lims = trialscan_helper(data_dir, pwrname, condname, channel);
% pwrlims = lims;
% 
% % 23 only
% condname = 'F023A159_F200A000';
% 
% figure(1)
% subplot(2, 2, 3)
% lims = trialscan_helper(data_dir, erpname, condname, channel);
% erplim = [erplim; lims];
% 
% figure(2)
% subplot(2, 2, 3)
% trialscan_helper(data_dir, pwrname, condname, channel)
% 
% % max both
% condname = 'F023A159_F200A016';
% 
% figure(1)
% subplot(2, 2, 4)
% trialscan_helper(data_dir, erpname, condname, channel)
% 
% figure(2)
% subplot(2, 2, 4)
% trialscan_helper(data_dir, pwrname, condname, channel)


%-----SUBFUNCTIONS BELOW-----
function data = trialscan_load(data_dir, filename, condname)
filename = [filename(1:23) condname filename(41:end)];
load([data_dir filename(42:end) '/' filename '.mat'])

function [xax, twod] = trialscan_data(data, channel)
twod = data.trial(channel, :, :);
twod = permute(twod, [3, 2, 1]);
try
    xax = data.time{1};
catch
    xax = data.freq{1};
end

% function clims = trialscan_helper(data_dir, filename, condname, channel)
% % load data
% filename = [filename(1:23) condname filename(41:end)];
% load([data_dir filename(42:end) '/' filename '.mat'])
% 
% % extract plotting data
% plotdat = data.trial(channel, :, :);
% plotdat = permute(plotdat, [3, 2, 1]);
% try
%     xax = data.time{1};
% catch
%     xax = data.freq{1};
% end
% 
% clims = find_lims(plotdat);
% 
% % plot
% imagesc(xax, [], plotdat, clims)


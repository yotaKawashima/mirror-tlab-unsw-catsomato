% This m-file was written Jul 2016 in order to investigate the reason that
% fieldtrip had such a field day with the timestamp. It was found that the
% data had uneven sampling, but we couldn't tell why. 
%
% So the purpose of this file is to plot the reported time resolution at
% the same time as the actual difference between timepoints.





data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/raw/';
name = 'Cat20110808_UtahRight-4_LFPStimSegs.mat'; % use the original file
loadme = [data_dir, name];

outname = 'C20110808_R04'; % need to change this line so that the names match up

clf

load(loadme, 'timestamp', 'tres')

ts1 = diff(timestamp);
plot(ts1, '.')

ts1_l = numel(ts1);

hold on; plot(1:ts1_l, tres*ones(1, ts1_l)); hold off


xlabel('sample number')
ylabel('difference from previous sample (s)')
legend({'data', 'expected'}, 'Orientation', 'horizontal', 'Location', 'SouthOutside')
title(outname, 'Interpret', 'none')

shg

print(gcf, '-dpng', [outname '_originaltimestampdataloss'])
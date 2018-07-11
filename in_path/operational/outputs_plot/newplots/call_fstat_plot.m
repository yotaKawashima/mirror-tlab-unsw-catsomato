data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/';
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';
filename = 'C20110808_R03';

dir_sec = [data_dir data_type];

options.load = false;
options.plotchannelscan = true;
options.plotchannel = {[]; [101]};
options.printtype = {'-dpng'};

%mvars = [];
%f_stats = [];

[mvars, f_stats] = fstat_plot(dir_sec, filename, options, mvars, f_stats);
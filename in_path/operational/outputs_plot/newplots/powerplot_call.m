% file to call the plot

options.frange = [0 250];
options.area = 'S1';
options.channel = 256;
options.version = 2016;
options.print = {'-dpng', '-depsc'};

options.condition = {'F023A159_F200A016'};
options.frequency = 23;

data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr';
filename = 'C20110808_R03';
%data = [];
%flag = 1;



[data, ~, ~] = powerplot(data_dir, filename, options, data);





data_dir = '/media/rannee/UNSW_Cat_Somatos/data/collated_data';

imgformat = '-depsc';


% power
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_avtrcon/';
collectallcats(data_dir, data_type, 'plotmean', 1)
print(1, imgformat, 'Figure2A')

% SNR
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr_avtrcon/';
collectallcats(data_dir, data_type, 'plotmean', 2, 'plotttest', 4)
print(2, imgformat, 'Figure2B')
print(4, imgformat, 'Figure2D')

% EP
data_type = 'epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_avtrcon/';
collectallcats(data_dir, data_type, 'plotmean', 3, 'plotttest', 5)
print(3, imgformat, 'Figure2C')
print(5, imgformat, 'Figure2E')
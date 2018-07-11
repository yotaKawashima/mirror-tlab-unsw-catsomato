area = 'S2';
chnum = 166;
erptype = 'epoched_dwnsmpl_biprref_evkresp';
pwrtype = 'epoched_dwnsmpl_biprref_evkresp_pwrspec_f0to250';

erpname = ['C20110808_R03_TStim_' area '_F023Axxx_F200Axxx_' erptype];
pwrname = ['C20110808_R03_TStim_' area '_F023Axxx_F200Axxx_' pwrtype];

data_dir = '/Users/ranneeli/Documents/MATLAB/ftformat_pipe/data/';

plot_trialscan(data_dir, erpname, pwrname, chnum)

print(1, '-dpng', ['trialscan_erp_' area '_ch' num2str(chnum, '%03i') '_' erpname])
print(3, '-dpng', ['trialscan_erp_' area '_ch' num2str(chnum, '%03i') '_' erpname '_colbar'])

print(2, '-dpng', ['trialscan_pwr_' area '_ch' num2str(chnum, '%03i') '_' pwrname])
print(4, '-dpng', ['trialscan_pwr_' area '_ch' num2str(chnum, '%03i') '_' pwrname '_colbar'])
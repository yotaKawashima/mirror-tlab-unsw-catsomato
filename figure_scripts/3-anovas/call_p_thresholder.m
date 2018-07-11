
loadme = '/media/rannee/UNSW_Cat_Somatos/data/included_datasets/C20110510_R05_S1_F023A1F000to250_P2_anoved_rsampsl_biprref_evkresp_cmtspwr_adatout.mat';

data_dir = '/media/rannee/UNSW_Cat_Somatos/data/collated_data/anoved_rsampsl_biprref_cmtspwr';
a = 1; 
cat_name = 'C20110510_R05';

filename_out = [cat_name '*S' num2str(a)];
q=0.05;

p_thresholder(data_dir, 'C20110808_R03', q, false)
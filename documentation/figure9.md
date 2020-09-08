# Figure 9 #

 Title: Nontagged frequencies between 50-150Hz in S2 are modulated by F2=200Hz vibratory amplitude.

 Show nontagged frequencies between 50Hz-150Hz were modulated by the vibration amplitude of F2=200Hz in only S2.

 Fig (A) --- F1 main effect f-statistis from ANOVA performed on logP in S1 and S2. <br />
 Fig (B) --- F2 main effect f-statistis from ANOVA performed on logP in S1 and S2.Not in matlab. <br />
 Fig (C) --- Interaction f-statistis from ANOVA performed on logP in S1 and S2.Time domain representatation of the bipolar-re-referenced LFP signal. <br />
 Fig (D) --- Exemplar nontagged VELogP responses between 50Hz-150Hz at bipolar channel 102 in S2 from Session 2-1. <br />
 Fig (E) --- Spatial mapping of HGP in S2. 
 

## File: plot_fstat_allfreq_yota_v4.m ##

**File discription** <br />
 Plot Fig(A) - (C) i.e. f-statistics from ANOVA performed on logP in S1 and S2. Note that ANOVA results of logP are the same as ones of VELogP.  

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly/

**Notes** 
+ Some channels in S1 gave invalid ANOVA results. This would be because of devision by 0 at logSNR computation (when responsese at neighboring frequenices are too small and matlab round down to 0). Need further confirmaton.

## File: S2_hgp_powerbycond_yota.m ##

**File discription** <br />
 Plot Fig (D) i.e. exemplar nontagged VELogP responses between 50Hz-150Hz at bipolar channel 102 in S2 from Session 2-1. Show mean and standard deviation accross 15 trials for each condition.

**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/

**Notes** 
+ []


## File: plot_S2_hgpinspace_yota.m ##

**File discription** <br />
  Plot Fig (E) i.e. spatial mapping of HGP in S2 from Session 2-1. Color encodes mean nontagged HGP across 15 trials for each condition.

**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/

**Notes** 
+ []



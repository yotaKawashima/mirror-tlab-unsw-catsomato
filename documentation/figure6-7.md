# Figure 6 #

**Description** 
 Title: Nontagged frequencies between 50-150Hz in S2 are modulated by F2=200Hz vibratory amplitude.

 Show nontagged frequencies between 50Hz-150Hz were modulated by the vibration amplitude of F2=200Hz in only S2.

 (A) --- F1 main effect f-statistis from ANOVA performed on logP in S1 and S2.
 (B) --- F2 main effect f-statistis from ANOVA performed on logP in S1 and S2.Not in matlab.
 (C) --- Interaction f-statistis from ANOVA performed on logP in S1 and S2.Time domain representatation of the bipolar-re-referenced LFP signal.


## File: plot_fstat_allfreq_yota_v4.m ##

**File discription** 
 Plot Fig 6A-C i.e. f-statistics from ANOVA performed on logP in S1 and S2. Note that ANOVA results of logP are the same as ones of VELogP.  

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly/

**Notes** 
+ Some channels in S1 gave invalid ANOVA results. Thse channels correspond to ones that data were not recorded at. (You can check it by running check_invalid_channels_v3.m.)


 
# Figure 7 #

 **Description** 
 Title: Exemplar nontagged frequencies responses between 50-150Hz modulated by F2=200Hz vibratory amplitude.

 (A) --- Exemplar nontagged VELogP responses between 50Hz-150Hz at bipolar channel 102 in S2 from Session 2-2.
 (B) --- Spatial mapping of HGP in S2. 
 

## File: S2_hgp_powerbycond_yota.m ##

**File discription** 
 Plot Fig 7A i.e. exemplar nontagged VELogP responses between 50Hz-150Hz at bipolar channel 102 in S2 from Session 2-2. Show mean and standard deviation accross 15 trials for each condition.

**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/

**Notes** 
+ []


## File: plot_S2_hgpinspace_yota.m ##

**File discription** 
  Plot Fig 7B i.e. spatial mapping of HGP in S2 from Session 2-2. Color encodes mean nontagged HGP across 15 trials for each condition.

**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp/

**Notes** 
+ []



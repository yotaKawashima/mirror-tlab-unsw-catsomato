# Figure 3-8 #

**Description** 
 Fig 3
 Title: Exemplar logSNR at f1=23Hz depend on the vibration amplitude of F1=23Hz. 

 Show mean and standard deviation across 15 trials for each condition.

 Fig (A) --- Exemplar responses at f1=23Hz at bipolar channel 131 in S1 from Session 2-1. 
 Fig (B) --- Spatial mapping of logSNR at f1=23Hz across all bipolar channels in S1.
 

 Fig 4
 Title: Exemplar logSNR at f2=200Hz depend on the vibration amplitude of F2=200Hz.

 Fig (A) --- Exemplar responses at f2=200Hz at bipolar channel 36 in S1 from Session 2-1.
 Fig (B) --- Spatial mapping of logSNR at f2=200Hz across all bipolar channels in S1.


 Fig 5
 Title: Exemplar logSNR at f1=23Hz depend on the vibration amplitude of F1=23Hz, F2=200Hz, and their interaction. 

 Fig (A) --- Exemplar responses at f1=23Hz at bipolar channel 158 in S1 from Session 2-1.
 Fig (B) --- Spatial mapping of logSNR at f1=23Hz across all bipolar channels in S1.


 Fig 6
 Title: Exemplar VELogP at f1=23Hz depend on the vibration amplitude of F1=23Hz.

 Fig (A) --- Exemplar responses at f1=23Hz at bipolar channel 176 in S1 from Session 2-1.
 Fig (B) --- Spatial mapping of VELogP at f1=23Hz across all bipolar channels in S1.


 Fig 7
 Title: Exemplar VELogP at f2=200Hz depend on the vibration amplitude of F2=200Hz.

 Fig (A) --- Exemplar responses at f2=200Hz at bipolar channel 122 in S1 from Session 2-1.
 Fig (B) --- Spatial mapping of VELogP at f2=200Hz across all bipolar channels in S1.


 Fig 8
 Title: Exemplar VELogP at f1=23Hz depend on the vibration amplitude of F1=23Hz, F2=200Hz, and their interaction.

 Fig (A) --- Exemplar responses at f1=23Hz at bipolar channel 43 in S1 from Session 2-1.
 Fig (B) --- Spatial mapping of VELogP at f1=23Hz across all bipolar channels in S1.


## File: get_pval_for_significant_channels.m ##

**File discription** 
 Print channels whose p values were smaller than FDR corrected threshold. 

**Data dependencies**
+ data/collated_data/anova_props/C20110808_R03\*snrsurr_adatain_adatout_pthresh.mat
+ data/collated_data/anova_props/C20110808_R03\*evkdpwr_adatain_adatout_pthresh.mat

**Notes** 
+ This is not necessary for plotting Fig3-8. 


## File: plot_examplepairs_yota_v2.m ##

**File discription** 
 Plot mean and standard deviation of responsese across 15 trials at a bipolar channel for each condidiotn. Also plot spatial map of mean responses at a frequency for each condition. 

**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr


# Figure 2 #

**Description** 
 Title: Results of the statistical analysis (two-way ANOVA on logP, logSNR and VELogP).

 Show % of channels that were deemed as significant according to two-way ANOVA for each freqency.

 Fig (A) --- % of channels showing siginficance on input amplitude of F1=23Hz main effect only.
 Fig (B) --- % of channels showing significance on input amplitude of F2=200Hz main effect only. 
 Fig (C) --- % of channels showing significance on the F1 and F2 main effects and F1-F2 interaction.
 

## File: call_p_thresholder_filenames.m ##

**File discription** 
 Call a function called p_thresholder_filenames.m, which calculates p value threshold with eeglab_fdr.m, and save the results to data/collated_data/anova_props/ directory.

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr

**Notes** 
+ Use eeglab_fdr.m.


## File: plot_anova_proportions_yota_v4.m ##

**File discription** 
 Compute and plot proprtions of channels showing siginificance on 1) only F1 main effect, 2) only F2 main effect, and 3) F1 and F2 main effects and their interaction.

**Data dependencies**
+ data/collated_data/anova_props/

**Notes** 
+ Run call_p_thresholder_filenames.m before runnig this file.
+ Some channels in S1 gave invalid ANOVA results. This would be because of devision by 0 at logSNR computation (when responsese at neighboring frequenices are too small and matlab round down to 0). Need further confirmaton.
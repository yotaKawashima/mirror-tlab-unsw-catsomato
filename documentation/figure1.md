# Figure 1 #

**Description** 
 Title: Neural recording, vibratory stimulation and our analysis scheme. 

 For each process, show exemplar responses at bipolar channel 156 in S1 from Session 2-2 (C20110808_R03). Specifically, plot mean and standard deviation across 15 trials from the max amplitude condition in the session i.e. [F1, F2]=[159, 16].

 (A) --- Not in matlab.
 (B) --- Not in matlab.
 (C) --- Time domain representatation of the bipolar-re-referenced LFP signal.
 (D) --- Power spectrum of the LFP signal.
 (E) --- Frequency domain representation of the LFP as logSNR. 
 (F) --- Frequency domain representation of the LFP as vibration evoked logPower (VELogP).

 
## File: summaryfigure.m ##

**File discription** 
 Load data first. Then, compute mean and standard deviation. Finally, plot them.
 
**Data dependencies**
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref 
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr
+ data/included_datasets/C20110808_R03/epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr

**Notes** 
+ t=0 is the start of the vibration.
+ Run in_path/operational/anova_analysis/pip_full_allcat.m in order to get data to plot beforehand.
+ Need to set directories in the heading.

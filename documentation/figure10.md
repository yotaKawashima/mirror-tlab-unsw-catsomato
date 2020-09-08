# Figure 10 #

 Title: Time course of the mean nontagged HGP in S2.

 Show that broadband (50Hz-150Hz) HGP are sustained after the onset transient. 

 Fig (A) --- Mean normalised power for HGP and frequencies of interest around HGB from the top 10% channels in S2. <br />
 Fig (B) --- Mean and Standard Error of the mean for HGP and F1=23Hz harmonics and Intermodulation around HGB from the top 10% channels in S2. <br />
 

## File: plottimecourse_all_yota_v3.m ##

**File discription** <br />
 Plot Fig (A) and Fig (B) i.e. mean normalised power for HGP and frequencies of interest around HGB from the top 10% channels in S2. This file performs the normalisation first. Specifically, we defined a baseline for each trial. As the baseline, we used mean power within a prestimulus period of -0.5s to -0.25s. (Here, HBW=2Hz corresponds to time window=0.5s.) After performing the normalisation per trial, compute and plot mean and standard deviation.

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_adatain_adatout_fstonly/
+ data/included_datasets/\*/epoched_rsampsl_biprref_cmtsgrm/

**Notes** 
+ When time course data are not computed previously, this file will automatically compute them
+ For Fig (B), you may want to change alpha values for shade. 


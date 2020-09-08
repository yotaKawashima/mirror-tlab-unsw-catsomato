# Figure 11 #

 Title: Waveforms of modeled processes.

 Fig (A) --- Waveforms in the time domain. <br />
 Fig (B) --- Waveforms in the frequency domain. <br />


## File: simulation_v2.m ##

**File discription** <br />
 Plot Fig (A) and (B). 

**Data dependencies**
+ NA

**Notes** 
+ []



# Figure 12 #

 Title: Exemplar model responses with logSNR observed at a channel.  

 Show exemplar model fitting results when fitting models to bipolar channel 43 in S1 from Session 2-1.

 Fig (A) --- Responses from four rectirication models. <br />
 Fig (B) --- Difference between rectification model logSNR and the observed logSNR. <br />
 Fig (C) --- Responses from two square wave models. <br />
 Fig (D) --- Difference between square wave model logSNR and the observed logSNR. <br />


## File: find_harmonics_IMs_channels_v3.m ##

**File discription** <br />
 Get top 10% channels to which we fit models. 

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr/

**Notes** 
+ Run this script before running other scripts related to model fitting.
+ This file saves 'Harmonics_IMs_top10.mat' in the script directory, the mat file that is required for model fitting analysis.


## File: find_harmonics_IMs_channels_v3.m ##

**File discription** <br />
 Get top 10% channels to which we fit models. 

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr/

**Notes** 
+ Run this script before running other scripts related to model fitting.
+ This file saves 'Harmonics_IMs_top10.mat' in the script directory, the mat file that is required for model fitting analysis.



## File: run_rect_batch.sh and run_sq_batch.sh ##

**File discription** <br />
  Slurm script to run matlab functions to fit models to top 10% channels

**Data dependencies**
+ figure_scripts/11-14-simulation/Harmonics_IMs_top10.mat

**Notes** 
+ Use this slurm file if you have access to a supercomputer with slurm system.
+ Replace ???? with your account, user name, and email address in the file and rect_top10_batch.sh and sq_top10_batch.sh.
+ Use run_rect_sq_with_matlab.m instead if you do not have access to a supercomputer. I have not tried this script, so I am not sure how long it would take. 
+ Run concatenate_fitting_results_from_massive.m after getting the results.

## File: concatenate_fitting_results_from_massive.m ##

**File discription** <br />
  Combine mat file containing model fitting results for each channel into one matlab structure. This script returns Resluts_Rect_top10_v2.m and Results_Sq_top10_v2.m, which are required for further analysis.

**Data dependencies**
+ figure_scripts/11-14-simulation/rect_top10_results_best_setting_v2/
+ figure_scripts/11-14-simulation/sq_top10_results_best_setting_v2/

**Notes** 
+ []  


## File: plot_exemplar_responses_with_model_v4.m ##

**File discription** <br />
 Plot Fig (A)-(D). 

**Data dependencies**
+ figure_scripts/11-14-simulation/Harmonics_IMs_top10.mat
+ figure_scripts/11-14-simulation/Results_Rect_top10_v2.mat

**Notes** 
+ Run this script after getting model fitting results.



# Figure 13 #

 Title: Comparison of fitting performances across rectification models and across square wave models.  

 Plot cumulative probability distribution of the minimum difference for comparision.

 Fig (A) --- Compare rectificatoin models in S1 and S2 separately. <br />
 Fig (B) --- Compare S1 and S2 for each rectification model. <br />
 Fig (C) --- Compare square wave models in S1 and S2 separately. <br />
 Fig (D) --- Compare S1 and S2 for each square model. <br />


## File: plot_model_results_v5.m ##

**File discription** <br />
 Plot Fig (A)-(D) with Kologorov-Smirnov test for Fig (B) and (D). 

**Data dependencies**
+ figure_scripts/11-14-simulation/Results_Rect_top10_v2.mat
+ figure_scripts/11-14-simulation/Results_Sq_top10_v2.mat

**Notes** 
+ []



# Figure 14 #

 Title: Akaike’s Information Criterion (AIC) for model comparison.  

 Compute AIC and plot for each somatosensory area.


## File: plot_AIC_v2.m ##

**File discription** <br />
 Compute AIC and plot Fig 14. 

**Data dependencies**
+ figure_scripts/11-14-simulation/Results_Rect_top10_v2.mat
+ figure_scripts/11-14-simulation/Results_Sq_top10_v2.mat

**Notes** 
+ []


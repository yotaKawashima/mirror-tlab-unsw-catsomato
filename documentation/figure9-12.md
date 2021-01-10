# Figure 9 #

**Description** 
 Title: Waveforms of modeled processes.

 (A) --- Waveforms in the time domain.
 (B) --- Waveforms in the frequency domain.


## File: simulation_v3.m ##

**File discription** 
 Plot Fig 9A and 9B. 

**Data dependencies**
+ NA

**Notes** 
+ []



# Figure 10 #

**Description** 
 Title: Exemplar channel’s observed and modelled responses.  

 Show exemplar model fitting results when fitting models to bipolar channel 43 in S1 from Session 2-2.

 (A) --- Responses from four rectirication models.
 (B) --- Difference between rectification model logSNR and the observed logSNR.


## File: find_harmonics_IMs_channels_v3.m ##

**File discription** 
 Get top 10% channels to which we fit models. 

**Data dependencies**
+ data/included_datasets/\*/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr/

**Notes** 
+ Run this script before running other scripts related to model fitting.
+ This file saves 'Harmonics_IMs_top10.mat' in the script directory, the mat file that is required for model fitting analysis.


## File: run_rect_batch.sh and run_hsq_batch.sh ##

**File discription** 
  Slurm script to run matlab functions to fit models to top 10% channels

**Data dependencies**
+ figure_scripts/9-12-simulation/Harmonics_IMs_top10.mat

**Notes** 
+ Use this slurm file if you have access to a supercomputer with slurm system.
+ Replace ???? with your account, user name, and email address in the file and rect_top10_batch.sh and hsq_top10_batch.sh.
+ Use run_rect_hsq_with_matlab.m instead if you do not have access to a supercomputer. I have not tried this script, so I am not sure how long it would take. 
+ Run concatenate_fitting_results_from_massive.m after getting the results.


## File: concatenate_fitting_results_from_massive.m ##

**File discription** 
  Combine mat file containing model fitting results for each channel into one matlab structure. This script returns Resluts_Rect_top10_v2.m and Results_HSq_top10_v2.m, which are required for further analysis.

**Data dependencies**
+ figure_scripts/9-12-simulation/rect_top10_results_best_setting_v2/
+ figure_scripts/9-12-simulation/hsq_top10_results_best_setting_v2/

**Notes** 
+ []  


## File: plot_exemplar_responses_with_model_v5.m ##

**File discription** 
 Plot Fig (A)-(B). 

**Data dependencies**
+ figure_scripts/9-12-simulation/Harmonics_IMs_top10.mat
+ figure_scripts/9-12-simulation/Results_Rect_top10_v2.mat

**Notes** 
+ Run this script after getting model fitting results.



# Figure 11 #

**Description** 
 Title: Comparison of fitting performances across rectification models.  

 Plot cumulative probability distribution of the minimum difference for comparision.

 (A) --- Compare rectificatoin models in S1 and S2 separately.
 (B) --- Compare S1 and S2 for each rectification model.


## File: plot_model_results_v6.m ##

**File discription** 
 Plot Fig 11A-11B with Kologorov-Smirnov test. 

**Data dependencies**
+ figure_scripts/9-12-simulation/Results_Rect_top10_v2.mat


**Notes** 
+ []



# Figure 12 #

**Description** 
 Title: Akaike’s Information Criterion (AIC) for 6 models for S1 and S2.  

 Compute AIC and plot for each somatosensory area.


## File: plot_AIC_v3.m ##

**File discription** 
 Compute AIC and plot Fig 12. 

**Data dependencies**
+ figure_scripts/9-12-simulation/Results_Rect_top10_v2.mat
+ figure_scripts/9-12-simulation/Results_HSq_top10_v2.mat
**Notes** 
+ []



# Sup Figure 7 #

**Description** 
 Title: Exemplar channel’s observed and modelled responses for half squaring models.  

 Show exemplar model fitting results when fitting models to bipolar channel 43 in S1 from Session 2-2.

 (A) --- Responses from two half squaring models.
 (B) --- Difference between half squaring model logSNR and the observed logSNR.


## File: plot_exemplar_responses_with_model_v5.m ##

**File discription** 
 Plot Sup Fig 4A-4B. 

**Data dependencies**
+ figure_scripts/9-12-simulation/Harmonics_IMs_top10.mat
+ figure_scripts/9-12-simulation/Results_HSq_top10_v2.mat

**Notes** 
+ Run this script after getting model fitting results.



# Sup Figure 8 #

**Description** 
 Title: Distribution of each coefficient in rectification models.

 Plot cumulative probability distribution of parameters for each model.

 (A) --- Results from S1 channels.
 (B) --- Results from S2 channels.
 

## File: plot_histogram_parameters_v3.m ##

**File discription** 
 Plot Sup Fig 6A-6B. 

**Data dependencies**
+ figure_scripts/9-12-simulation/Results_Rect_top10_v2.mat

**Notes** 
+ []

# Sup Figure 9 #

**Description** 
 Title: Comparison of fitting performances across half squaring models.  

 Plot cumulative probability distribution of the minimum difference for comparision.

 (A) --- Compare half squaring models in S1 and S2 separately.
 (B) --- Compare S1 and S2 for each half squaring.


## File: plot_model_results_v6.m ##

**File discription** 
 Plot Sup Fig 5A-5B with Kologorov-Smirnov test. 

**Data dependencies**
+ figure_scripts/9-12-simulation/Results_HSq_top10_v2.mat

**Notes** 
+ []


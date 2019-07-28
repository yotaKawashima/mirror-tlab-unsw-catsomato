# Figure 3 #

**Description:** Results of 2-way ANOVA at FOIs

## File: call_p_thresholder.m ##

**File description:** What it says on the tin.

**Data dependencies:**

+ anoved_rsampsl_biprref_evkresp_cmtspwr_adatout for all cats

**Notes:** 

+ Output data type: anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_pthresh
+ Need to set directories in the heading before running.

**TODO:** 

+ Move to an appropriate directory
+ []
+ Add script dependencies to wiki
+ Automatic adding of the git repo to path at matlab startup


## File: plot_anova_proportions ##

**File description:** Plots figure 3

**Data dependencies:**

+ anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_pthresh for all cats (created from anoved_rsampsl_biprref_evkresp_cmtspwr_adatout using `call_p_thresholder`)

**Notes:** 

+ Produces two alternate versions of figure 3B, one of which only shows frequencies of interest, and the other shows all frequencies.
+ Last version before version control is /media/rannee/UNSW_Cat_Somatos/scripts/Nov17/fig3_20171113_b.m

**TODO:** 

+ NA
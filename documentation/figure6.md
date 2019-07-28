# Figure 6 #

**Description:** Non-tagged high gamma power

## File: plot_fstat_allfreq.m ##

**File description:** plots A and B

**Data dependencies:**

+ anoved_rsampsl_biprref_evkresp_cmtspwr_adatout_fstonly (contains data from all cats)

**Notes:** 

+ Need to set directories in the heading before running.
+ Last version before version control is /media/rannee/UNSW_Cat_Somatos/scripts/Nov17/fstat_fig6_plotter_20171126.m
    + That script was based on /media/rannee/UNSW_Cat_Somatos/scripts/28Jul17/fstat_fig4_plotter.m

**TODO:** 

+ Automatic adding of the git repo to path at matlab startup
+ Figure filenames
+ Work out higher level data dependencies
+ Add script dependencies to wiki

## File: S2_hgp_powerbycond.m ##

**File description:** Plots C

**Data dependencies:**

+ epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp for C20110808_R03 only

**Notes:** 

+ Input is bi-polar channel number
+ Need to set directories in the heading before running.
+ Last version before version control is /media/rannee/UNSW_Cat_Somatos/scripts/Nov17/figure3_hgp_14Nov.m

**TODO:** 

+ Automatic adding of the git repo to path at matlab startup
+ Add script dependencies to wiki


## File: plot_S2_hgpinspace.m ##

**File description:** plots D

**Data dependencies:**

+ epoched_rsampsl_biprref_evkresp_cmtspwr_evkdpwr_hgpcomp for C20110808_R03 only

**Notes:** 

+ Colour bar is truncated and upside down.
+ []
+ Need to set directories in the heading before running.
+ Last version before version control is /media/rannee/UNSW_Cat_Somatos/scripts/Nov17/figure7b_20171114_a.m
    + This plotted both S1 and S2.

**TODO:** 

+ Add the green square that corresponds to the channel in C
+ Add the lower level data dependencies
+ []
+ Automatic adding of the git repo to path at matlab startup
+ Add script dependencies to wiki
# Code Description

The following R script files with are
used to run an experiment and to save the results. 
The below lists the main files in the order they should be run. 

### Main code files

* `01_deconvolve_daily_reinfects.R` and `01_deconvolve_weekly_reinfects.R` are for deconvolving daily, 
  weekly or biweekly reinfections back to the approximate date of infection onset.
* `02_ready_sero_for_ssmod.R` readies the commercial lab and blood donor seroprevalence survey data for 
  use in the state space model.
* `03_state_space_model.R` executes the leaky immunity state space model to obtain
  the inverse reporting ratios.
* `04_gp_inv_ratios.R` runs the Gaussian process code to generate multiple sets of plausible inverse reporting ratios
  and infection trajectories per state.
* `05_estimate_infections_from_ww.R` is the main script to estimate the variant-specific shedding rates and inverse 
  reporting ratios over Omicron.
* `06_estimate_infections_from_ww_sensitivity.R` runs the sensitivity analysis of the infection estimates to the 
  assumed shedding profile.
* `07_ww_growth_rates.R` generates the time-varying growth rates per state.
* `08_ww_rt_estimation.R` generates the time-varying instantaneous reproduction numbers per state.
* `09_ww_rt_sensitivity.R` evaluates the sensitivity of Rt to varying generation times and serial intervals.
* `10_ww_growth_rates_sensitivity.R` evaluates the sensitivity of growth rates to varying window sizes.

### Supplementary files 
The following supplementary files for the deconvolution-based experiments are available:
 * `supporting-files-for-deconvolutions` holds several R and C++ functions
  that are necessary to load before attempting the deconvolutions by state. 

### Generate results

The generate-numeric-results fodler contains an .Rmd file to produce the numerical results in the paper. 
Note that the data used to generate these results is 
contained in the data folder, which also contains the public data inputs and the derived 
state-specific `convolution-mat-list.rds` files that are central to the deconvolutions.


## Credit
This repository follows and builds off of the code and other materials from
[stat-sci-nowcast](https://github.com/cmu-delphi/stat-sci-nowcast/) and 
[latent-infections](https://github.com/cmu-delphi/latent-infections/).

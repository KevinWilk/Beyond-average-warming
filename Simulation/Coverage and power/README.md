## Overview 

#### Coverage_and_power_of_Test.R
- Runs a simulation data setup 
- Estimates quantile via the dependent/independent Multiplier Bootstrap with N^* = 1000
- constructs confidence bands
- after 1000 repetition: calculates empirical coverage rate
- after 1000 repetition: calculates emirical power (under Alternative)
  
#### Run_simluation.R
- runs ‘Coverage_and_power_of_Test.R’ for different parameters
  
#### Run_simulation_60.sh
- runs ‘Run_simluation.R‘ with 60 cores (MaRC3a, parallelization)
#### Run_simulation_120.sh
- runs ‘Run_simluation.R‘ with 120 cores (MaRC3a, parallelization)
  
#### Folder:
##### Under H0: p = 25, p = 50, p = 75
- saved results of empirical quantile of 1000 repititions (n_rep_1000_emp(1-alpha).rds), empirical coverage rate (dep/ind_n_rep_1000_emp(1-alpha).rds), list of 1000 estimated quantiles (dep/ind_n_rep_1000_q(1-alpha)_list.rds)
##### Under H1: only for p = 25
- saved results of empirical quantile of 1000 repititions (n_rep_1000_emp(1-alpha).rds), empirical coverage rate (dep/ind_n_rep_1000_emp(1-alpha).rds), empirical power (dep/ind_n_rep_1000_power(1-alpha).rds), list of 1000 estimated quantiles (dep/ind_n_rep_1000_q(1-alpha)_list.rds)

## Overview 

#### `Coverage_and_power_of_Test.R`
- Runs a simulation data setup 
- Estimates quantile via the dependent/independent Multiplier Bootstrap with N^* = 1000
- constructs confidence bands
- after 1000 repetition: calculates empirical coverage rate
- after 1000 repetition: calculates emirical power (under Alternative)
  
#### `Run_simluation.R`
- runs ‘Coverage_and_power_of_Test.R’ for different parameters

   
##### Under H0: p = 25, p = 50, p = 75 
- applied dependent Multiplier Bootstrap (`N* = 1000`) to estimate 90%-, 95%- and 99%-quantiles and saved as list of 1000 estimated quantiles (`dep_n_rep_1000_q_list.rds`)
- saved results of empirical quantiles of 1000 repititions (`n_rep_1000_emp.rds`) as list
- saved empirical coverage rates (`dep_n_rep_1000_coverage_list.rds`) as list
##### Under H1: p = 25, p = 50, p = 75
- applied dependent Multiplier Bootstrap (`N* = 1000`) to estimate 90%-, 95%- and 99%-quantiles and saved as list of 1000 estimated quantiles (`dep_n_rep_1000_q_list.rds`)
- saved results of empirical quantiles of 1000 repititions (`n_rep_1000_emp.rds`) as list
- saved empirical coverage rates (`dep_n_rep_1000_coverage_list.rds`) as list
- saved empirical power (`dep_n_rep_1000_power.rds`) as list

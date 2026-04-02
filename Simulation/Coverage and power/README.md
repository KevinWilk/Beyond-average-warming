## Overview 

#### Coverage_and_power_of_Test.R
- Runs a simulation data setup 
- Estimates quantile via the dependent/independent Multiplier Bootstrap with N^* = 1000
- constructs confidence bands
- after 1000 repitition: calculates empirical coverage rate
- after 1000 repitition: calculates emirical power (Under Alternative)
  
#### Run_simluation.R
- runs Coverage_and_power_of_Test.R for different parameters
  
#### Run_simulation_60.sh
- 
#### Run_simulation_120.sh
- 
#### Folder:
##### Under H0: p 25, p 50, p 75
- saved results of : empirical coverage rate (dep/ind_n_rep_1000_emp(1-alpha).rds), empirical quantile of 1000 repititions, list of 1000 estimated quantiles
##### Under H1
- saves empirical coverage rate

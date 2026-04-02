## Overview 

#### `LagCov_bandwidth_selection.R`: Bandwidth selection for lagged covariance kernel
- `bw.list.p.25/50/75` and `bw.list.pd.100`: via hv-crossvalidation in a 5-fold framework <br>
  -> saved as `lag_Gamma_s_25/50/75_Bandwidth.rds` in *Results* <br>
  -> saved as `lag_Gamma_d_100_Bandwidth.rds` in *Results/Delta = 0* (under H0: delta constant) and *Results/Delta != 0* (under H0: delta not constant)
  
- `bw.sample.list.p.25/50/75` and `bw.sample.list.pd.100`: saved selected bandwidth of length 1000 <br> 
  -> saved as `lag_Gamma_s_25/50/75_sample.rds` in Results
  
- `bw.optim.list.p.25/50/75` and `bw.optim.list.pd.100`: via true underlying lagged covariance kernel <br>
  -> saved as `lag_Gamma_s_25/50/75_Bandwidth_optim.rds` in *Results* <br>
  -> saved as `lag_Gamma_d_100_Bandwidth_optim.rds` in *Results/Delta = 0* (under H0: delta constant) and *Results/Delta != 0* (under H0: delta not constant)
  
  
#### Mean_bandwidth_selection.R: Bandwidth selection for mean and difference function
- `delta.hv.25/50/75` and `dense.hv.100`: via hv-crossvalidation in a 5-fold framework
- `delta.optim.25/50/75` and `dense.optim.100`: via true underlying difference and mean function <br>
  -> saved as `delta_constant_TRUE_Bandwidth.rds` in *Results* (under H0: delta constant) <br>
  -> saved as `delta_constant_FALSE_Bandwidth.rds` in *Results* (under H1: delta not constant)
  
#### `LagCov_bandwidth_selection_cores_60.sh`
- runs `LagCov_bandwidth_selection.R` with 60 cores (MaRC3a, parallelization)
  
#### `Mean_bandwidth_selection_cores_60.sh`
- runs `Mean_bandwidth_selection.R` with 60 cores (MaRC3a, parallelization)
  


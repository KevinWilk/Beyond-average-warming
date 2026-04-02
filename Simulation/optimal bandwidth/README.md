## Overview 

#### LagCov_bandwidth_selection.R: Bandwidth selection for lagged covariance kernel
- bw.list.p.25/50/75 and bw.list.pd.100: via hv-crossvalidation in a 5-fold framework
- bw.sample.list.p.25/50/75 and bw.sample.list.pd.100: saved selected bandwidth of length 1000
- bw.optim.list.p.25/50/75 and bw.optim.list.pd.100: via true underlying lagged covariance kernel 
 
#### Mean_bandwidth_selection.R: Bandwidth selection for mean and difference function
- delta.hv.25/50/75 and dense.hv.100: via hv-crossvalidation in a 5-fold framework
- delta.optim.25/50/75 and dense.optim.100: via true underlying difference and mean function

#### LagCov_bandwidth_selection_cores_60.sh
- runs ‘LagCov_bandwidth_selection.R‘ with 60 cores (MaRC3a, parallelization)
  
#### Mean_bandwidth_selection_cores_60.sh
- runs ‘Mean_bandwidth_selection.R‘ with 60 cores (MaRC3a, parallelization)
  
#### Folder:
##### Results
-

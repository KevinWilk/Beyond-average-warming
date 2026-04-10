## Overview 

#### `function_bandwidth_selection.R` cotains

- `k.fold.hc.cv()`:     K-Fold hv-block cross validation for univariate local polynomial estimator
- `bw_month()`:         runs `k.fold.hc.cv()` on densely (2000-2025) observed data sets and sparsly (1952-1972) observed residuals for each season with `K=5`
- `cov.k.fold.hc.cv()`: K-Fold hv-block cross validation for bivariate local polynomial estimator
- `cov.bw.month()`:     runs `cov.k.fold.hc.cv()` for each season with `K=5`


#### `bandwidth_selection.R`
##### for mean/difference function
-  `Bandwidths`: via hv-crossvalidation in a K-fold framework for each season (K = length of time period)
-> saved as  `bw_city.rds` in _Results_ 
##### for covariance kernel (lag 0)
-  `cov.Bandwidths`: via hv-crossvalidation in a K-fold framework for each season (K = length of time period)
-> saved as  `bw_Gamma_s_city.rds` (sparse set) or `bw_Gamma_d_city.rds` (dense set) in _Results_ 

   

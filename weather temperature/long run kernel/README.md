## Overview 

#### `function_test.R` cotains
- `inference.autocovariance.test()`: Tests H0: Gamma( , ;b) = 0, where |b| = 1,...,12 
  
#### `lagged_covariance.R` cotains
- calculation of lagged covariance matrix for b = 0,...,12
  -> saved as `list_Gamma_s_city.rds` (sparse: 1952 - 1972) or `list_Gamma_d_city.rds` (dense: 2000 - 2025) in _Results_
- Test results of inference autocovariance for each b = 1,...,12
  -> saved as `test_max_lag_s_city.rds` (sparse: 1952 - 1972) or `test_max_lag_d_city.rds` (dense: 2000 - 2025) in _Results_
- Maximum lag of long run covariance kernel
  -> saved as `max.lag` in _Results_


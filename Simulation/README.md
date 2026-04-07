## Overview 

#### `Simulation_Mod_biLocPol.R`
- Grid setting `(1:j-0.5)/p` similar to Hajo Holzmann and Max Berger [arxiv](https://arxiv.org/abs/2407.13641)

#### `function_bandwidth_selection.R` cotains
- `h.optim`:          Optimal bandwidth selection with known function for univariate local polynomial estimator
- `h.optim.cov`:      Optimal bandwidth selection with known function for bivariate local polynomial estimator
- `k.fold.hc.cv`:     K-Fold hv-block cross validation for univariate local polynomial estimator 
- `k.fold.hc.cv.cov`: K-Fold hv-block cross validation for bivariate local polynomial estimator
  
#### `function_data_generating.R` cotains
- `mu`,`mu_d` (mean functions) and `delta` (difference function)
- `sim.d.OU`: discrete simulation of processes from section 4
- `cov.d.OU`: True lagged covariance function (default: lag = 0)
- `LR.cov.d.OU`: True long run covariance function

#### `function_evaluation.R` cotains
-  `est.results`: estimation of mean functions and differnce function
-  `q.MB`: prepares arguments (lists `ls` and `ld`,  `cov`,  `depend = TRUE/FALSE`,  `int= TRUE/FALSE`) for Multiplier Bootstrap function `MB` and runs `MB`-function `B=1000`-times
-  `MB`: indendent/dependent Multiplier Bootstrap

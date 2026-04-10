## Overview 

#### `functions.R` cotains
-  `est.results()`: estimation of mean functions (1952-1972 and 2000-2025) and difference function on 15-minutes grid (00:00 to 24:00) for each month 
- `cov.weights()`:  calculates weights of lagged bivariate local polynomial estimator for each month 
- `test.int()`:     significance test of averaged daily temperature (integral of delta) for each month
  
- `P.Cov()`: calculates linear projection (P) of lagged covariance matrix (_full_ or _diagonal_)
- `q.month()`: runs `q.MB()` for each month 
- `q.MB()`: prepares arguments (lists `ls` and `ld`,  `cov`,  `depend = TRUE/FALSE`,  `int= TRUE/FALSE`) for Multiplier Bootstrap function `MB` and runs `MB()`-function `B=1000`-times
-  `MB()`: indendent/dependent Multiplier Bootstrap
-  `CB()`: constructs simultaneous confidence bands for each month
  


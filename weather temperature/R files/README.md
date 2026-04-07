## Overview 

#### `function_bandwidth_selection.R` cotains
- `k.fold.hc.cv()`: K-Fold hv-block cross validation for univariate local polynomial estimator
- `bw_month()`: runs `k.fold.hc.cv()` on densly (2000-2025) sampled data sets and sparsly (1952-1972) sampled residuals
- ```latex
$$
Y^{[s]} - {\hat{Y}}^{[d]}
$$
-  `est.results()`: estimation of mean functions (1952-1972 and 2000-2025) and difference function on 15-minutes grid (00:00 to 24:00)
  
- `cov.k.fold.hc.cv()`: K-Fold hv-block cross validation for bivariate local polynomial estimator
- `cov.bw.month()`:
- `cov.weights()`:
  
- `test.int()`: for each month
  
- `P.Cov()`:
- `q.month()`:
- `q.MB()`: prepares arguments (lists `ls` and `ld`,  `cov`,  `depend = TRUE/FALSE`,  `int= TRUE/FALSE`) for Multiplier Bootstrap function `MB` and runs `MB()`-function `B=1000`-times
-  `MB()`: indendent/dependent Multiplier Bootstrap
-  `CB()`:
  



This repository contains all numerical results and code implementations associated with the paper:

# Beyond average warming: Two-sample inference for dense-sparse functional data reveals changes in intraday temperature patterns

<!-- badges: start -->
<!-- badges: end -->

Kevin Wilk and Hajo Holzmann [arxiv](http://arxiv.org/abs/2605.27335). 

Note the file `Mod_BiLocPol.R` is based on the package `biLocPol` of which you can get the development version from [GitHub](https://github.com/mbrgr/biLocPol) with:

``` r
# install.packages("devtools")
devtools::install_github("mbrgr/biLocPol")
# library(biLocPol)
```
Important: `Mod_biLocPol.R` contains apdapted and extended functions of this package.
``` r
source("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/Mod_biLocPol.R")
```

## Comments

Note that the calculation of the weights, the simulation and the application are paralallized with the “future.apply” and “future” package and performed using 
the Marburg Compute Cluster (MaRC3a) at the Philipps-Universität Marburg [Link](https://www.uni-marburg.de/en/hrz/services/high-performance-computing).

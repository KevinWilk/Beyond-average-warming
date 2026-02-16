# Beyond the positive drift

This repository contains all numerical results and code implementations associated with the paper:

# Beyond the positive drift: Comparing historical and current daily temperature patterns based on two sample statistics for unbalanced dense-sparse functional data

<!-- badges: start -->
<!-- badges: end -->

Kevin Wilk and Hajo Holzmann [arxiv]( ). 

Note the this code needs the package “biLocPol” of which you can get the development version from [GitHub](https://github.com/mbrgr/biLocPol) with:

``` r
# install.packages("devtools")
devtools::install_github("mbrgr/biLocPol")
# library(biLocPol)
```
Important: “Mod_biLocPol.R“ contains apdapted and extended functions of this package.


## Structure


## Comments

Some of the “.Rdata” files are large since they contain the actual weights of some of the estimations. Note that the calculation of the weights and the simulations are paralallized with the “future.apply” and “future” package and performed using 
the Marburg Compute Cluster (MaRC3a) at University of Marburg [Link](https://www.uni-marburg.de/en/hrz/services/high-performance-computing).

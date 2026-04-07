library(crayon)
library(ggplot2)
library(reshape2)
library(locpol)
library(biLocPol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(lubridate)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)


source("Mod_biLocPol.R")
source("weather temperature/functions.R")
source("weather temperature/bandwidth selection/function_bandwidth_selection.R")





# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 64 * 1024^3)                ##   
plan(multisession, workers = 20) # MaRC3a: partition=mqtest  ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################



########################################################################################################################
################ Covariance kernel: Cross Validation on sparsely observed data set (1952-1972) #########################
########################################################################################################################



for(k in 1:4){

  ############################################################################
  #        Change to: 1 (Berlin)                                            ##
  #                   2 (Frankfurt am Main)                                 ##
  #                   3 (Hamburg)                                           ##
  #                   4 (Munich)                                            ##
  data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData")) ##
  ############################################################################



  cov.Bandwidths = cov.bw.month(data.sparse.days.5h, h.seq = seq(0.1,0.4,0.01), period.n = 24)
  file = paste0("weather temperature/bandwidth selection/Results/bw_Gamma_s_",data.example[[k]],".rds")
  saveRDS(cov.Bandwidths,file)
  print("save done")

}
                          
plan(sequential)








########################################################################################################################
################ Covariance kernel: Cross Validation on densely observed data set (2000-2025) ##########################
########################################################################################################################


### Only 30 hours data set possilbe due to "Killed" for 34 hours dataset (60 cores).

# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 64 * 1024^3)                ##   
plan(multisession, workers = 20) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################


for(k in 1:4){
  
  ############################################################################
  #        Change to: 1 (Berlin)                                            ##
  #                   2 (Frankfurt am Main)                                 ##
  #                   3 (Hamburg)                                           ##
  #                   4 (Munich)                                            ##
  data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData")) ##  
  ############################################################################
  
  

  cov.Bandwidths = cov.bw.month(data.dense.days.3h, h.seq = seq(0.06,0.24,0.01), period.n = 144)
  file = paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[k]],".rds")
  saveRDS(cov.Bandwidths,file)
  print("save done")

}

plan(sequential)

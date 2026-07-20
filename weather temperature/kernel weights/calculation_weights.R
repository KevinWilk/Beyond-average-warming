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

########################################################################################################################
################ Lagged kernel: Calculating weights ####################################################################
########################################################################################################################





############################################################################
data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      ##
load(paste0("weather temperature/data sets/",data.example[[1]],".RData")) ##
############################################################################




# weights only depends on bandwidth for fixed design grid

cov.Bandwidths.1 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_s_",data.example[[1]],".rds"))
cov.Bandwidths.2 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_s_",data.example[[2]],".rds"))
cov.Bandwidths.3 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_s_",data.example[[3]],".rds"))
cov.Bandwidths.4 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_s_",data.example[[4]],".rds"))

bw.s = unique(c(t(cov.Bandwidths.1),t(cov.Bandwidths.2),t(cov.Bandwidths.3),t(cov.Bandwidths.4)))

cov.Bandwidths.1 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[1]],".rds"))
cov.Bandwidths.2 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[2]],".rds"))
cov.Bandwidths.3 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[3]],".rds"))
cov.Bandwidths.4 = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[4]],".rds"))

bw.d = unique(c(t(cov.Bandwidths.1),t(cov.Bandwidths.2),t(cov.Bandwidths.3),t(cov.Bandwidths.4)))




# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 60 * 1024^3)                ##   
plan(multisession, workers = 120) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################

p.eval = 137

# for 1952 - 1972 (sparse)

for(bw in bw.s){

  w  = local.polynomial.weights(35,  bw, p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "less", eval.type = "full", parallel.environment = F)
  saveRDS(w, file = "weather temperature/kernel weights/full_w_s_lag0.rds")
  rm(w); gc()
  
  print(paste0("done: lag 0"))
  
  for(l in 1:12){
    bw.lag = round(bw*1.1^l, digits = 2)
    w.lag  = local.polynomial.weights(35, bw.lag, p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)
    saveRDS(w.lag, file = paste0("weather temperature/kernel weights/full_w_s_lag",l,".rds"))
    
    print(paste0("done: lag ", l))
    
    rm(w.lag); gc()
  }
}



# for 2000 - 2025 (dense)

for(bw in bw.d){
  
  w  = local.polynomial.weights(205,  bw, p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "less", eval.type = "full", parallel.environment = F)
  
  if(bw < 0.1){ saveRDS(w, file = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%03d", as.integer(bw * 100)),".rds"))
  }else{        saveRDS(w, file = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%02d", as.integer(bw * 10)),".rds"))
  }
  rm(w); gc()
  
  print(paste0("done: lag 0"))
  
  for(l in 1:12){
    
    bw.lag = round(bw*(30/34)*1.1^l, digits = 2)  # transformation from 30h time period to 34h
    w.lag  = local.polynomial.weights(205, bw.lag, p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)
    
    if(bw < 0.1){ saveRDS(w.lag, file = paste0("weather temperature/kernel weights/full_w_d_lag",l,"_",sprintf("%03d", as.integer(bw * 100)),".rds"))
    }else{        saveRDS(w.lag, file = paste0("weather temperature/kernel weights/full_w_d_lag",l,"_",sprintf("%02d", as.integer(bw * 10)),".rds"))
    }
    
    print(paste0("done: lag ", l))
    
    rm(w.lag); gc()
    
  }
  
}


plan(sequential)















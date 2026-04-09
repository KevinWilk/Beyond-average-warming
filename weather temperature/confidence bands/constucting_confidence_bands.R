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
plan(multisession, workers = 20) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################






for(k in 1:4){


  
  
  ############################################################################
  k = 1  # Change to: 1 (Berlin)                                            ##
  #                   2 (Frankfurt am Main)                                 ##
  #                   3 (Hamburg)                                           ##
  #                   4 (Munich)                                            ##
  data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData")) ##
  ############################################################################



  # Loading estimated lag covariance kernels for lag = 0,...,12 
  cov.s.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_s_",data.example[[k]],".rds")) 
  cov.d.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_d_",data.example[[k]],".rds")) 

  # Loading maximum lag of long run kernel 
  max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds")) 

  # Loading bandwidths for mean and difference estimation
  Bandwidths = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_",data.example[[k]],".rds"))
  
  # Mean and difference function estimation
  est = est.results(data.s.34h,data.d.34h,Bandwidths,from = 1,to = 12)



  start.24 = 21;end.24 = 117

  #######################################
  # Calculating variance for each month #
  #######################################

  lr.var.d.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,
                                                      lapply(1:(max.lag[m,2]+1),  
                                                             function(j){if(j >= 2){
                                                               (cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                }else{cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
                                                    }))
  lr.var.s.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,
                                                      lapply(1:(max.lag[m,2]+1),  
                                                             function(j){if(j >= 2){
                                                               (cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                }else{cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
                                                    }))
  
  ################################################################
  # Integral of delta test with construction of confidence bands #
  ################################################################
  
  set.seed(2004)
  integral.test = test.int(data.s.34h, data.d.34h, unique(est$delta_int$ESTIMATE), lr.var.s.month, lr.var.d.month, from = 1, to = 12, test = "two-sided", alpha = 0.95)
  unique(est$delta_int$ESTIMATE)
  
  integral.conf       = tibble(TIME = rep(hms::as_hms(c(as.POSIXct("1970-01-01 00:00:00"),as.POSIXct("1970-01-01 23:59:59"))),times = 12))
  integral.conf$UP    = rep(unlist(integral.test$confInterval)[seq(2,24,2)], each = 2)
  integral.conf$LO    = rep(unlist(integral.test$confInterval)[seq(1,24,2)], each = 2)
  integral.conf$MONTH =  factor(rep(month.name , each = 2), level = month.name)

  file = paste0("weather temperature/confidence bands/Results/CB_integral_",data.example[[k]],".rds")
  saveRDS(integral.conf, file = file)


  
  
  
  
  
  #####################################################################
  # Calculating diagonal of long run covariance matrix for each month #
  #####################################################################

  lr.Gamma.d.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                            function(j){if(j >= 2){
                                                              (cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                }else{cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                        }))
                                              return(diag(lr.cov))})
  
  lr.Gamma.s.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                            function(j){if(j >= 2){
                                                              (cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                }else{cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                        }))
                                               return(diag(lr.cov))})

  lr.Gamma = lapply(1:12, function(m) { sample.dense  = data.d.24h.hourly %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d.24h.hourly)[2])
                                        sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
                                        sample.sparse = data.s.24h %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s.24h)[2])
                                        sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
                                        return(lr.Gamma.s.month[[m]] + dim(sample.sparse)[1] / dim(sample.dense)[1] * lr.Gamma.d.month[[m]])})

  #####################################################################################
  # Dependent Multiplier Bootstrap on non centered difference function for each month #
  #####################################################################################

  set.seed(2004)
  q.list = q.month(data.s.34h, data.d.34h, Bandwidths, est, lr.Gamma, from = 1, to = 12, alpha = 0.95, B = 10000, depend = T, int = F)
  file = paste0("weather temperature/confidence bands/Results/quantile_",data.example[[k]],".rds")
  saveRDS(q.list, file = file) 

  
  # Skip: 
  # q.list = readRDS(paste0("weather temperature/confidence bands/Results/quantile_",data.example[[k]],".rds"))

   
  #################################################################################
  # Difference function test (not centered) with construction of confidence bands #
  #################################################################################

  delta.conf = CB(data.s.34h, est, lr.Gamma, q.list[1,])
  file = paste0("weather temperature/confidence bands/Results/CB_delta_",data.example[[k]],".rds")
  saveRDS(delta.conf, file = file) 
  
  
  
  
  
  
  
  
  
  

  ##############################################################################################
  # Calculating diagonal of linear projection (P) of long run covariance matrix for each month #
  ##############################################################################################

  lr.P.Gamma.d.month = lapply(1:12, function(m) {
    
    start = 21; end = 116
    
    lr.cov.24 = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                  function(j){if(j >= 2){
                                    (cov.d.list[[m]][[j]][start:end,start:end]+t(cov.d.list[[m]][[j]][start:end,start:end]))*(1-(j-1)/(max.lag[m,2]+1))
                                      }else{cov.d.list[[m]][[j]][start:end,start:end]} 
                                              }))
    start = 21; end = 117
    
    lr.cov.25 = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                  function(j){if(j >= 2){
                                    (cov.d.list[[m]][[j]][start:end,start:end]+t(cov.d.list[[m]][[j]][start:end,start:end]))*(1-(j-1)/(max.lag[m,2]+1))
                                      }else{cov.d.list[[m]][[j]][start:end,start:end]} 
                                              }))
  
    return(diag(lr.cov.25) - 2*rowMeans(lr.cov.25) + mean(lr.cov.24)  )})


  lr.P.Gamma.s.month = lapply(1:12, function(m) {
    
    start = 21; end = 116
    
    lr.cov.24 = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                  function(j){if(j >= 2){
                                    (cov.s.list[[m]][[j]][start:end,start:end]+t(cov.s.list[[m]][[j]][start:end,start:end]))*(1-(j-1)/(max.lag[m,2]+1))
                                      }else{cov.s.list[[m]][[j]][start:end,start:end]} 
                                              }))
    start = 21; end = 117
    
    lr.cov.25 = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                  function(j){if(j >= 2){
                                    (cov.s.list[[m]][[j]][start:end,start:end]+t(cov.s.list[[m]][[j]][start:end,start:end]))*(1-(j-1)/(max.lag[m,2]+1))
                                      }else{cov.s.list[[m]][[j]][start:end,start:end]} 
                                              }))
    
    return(diag(lr.cov.25) - rowMeans(lr.cov.25) - colMeans(lr.cov.25) + mean(lr.cov.25)  )})


  P.Gamma.lr = lapply(1:12, function(m) { sample.dense  = data.d.24h.hourly %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d.24h.hourly)[2])
                                          sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
                                          sample.sparse = data.s.24h %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s.24h)[2])
                                          sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
                                          return(lr.P.Gamma.s.month[[m]] + dim(sample.sparse)[1] / dim(sample.dense)[1] * lr.P.Gamma.d.month[[m]])})


  #################################################################################
  # Dependent Multiplier Bootstrap on centered difference function for each month #
  #################################################################################

  set.seed(2004)
  q.P.list = q.month(data.s.34h, data.d.34h, Bandwidths, est, P.Gamma.lr, from = 1, to = 12, alpha = 0.95, B = 10000, depend = T, int = T) 
  file = paste0("weather temperature/confidence bands/Results/quantile_P_",data.example[[k]],".rds")
  saveRDS(q.P.list, file = file) 
  
  # Skip: 
  # q.P.list = readRDS(paste0("weather temperature/confidence bands/Results/quantile_P_",data.example[[k]],".rds"))
  
  
  
  ###########################################################################
  # Centered difference function test with construction of confidence bands #
  ###########################################################################
  
  centered.delta.conf = CB(data.s.34h, est, P.Gamma.lr, q.P.list[1,], center = T)
  file = paste0("weather temperature/confidence bands/Results/CB_centered_delta_",data.example[[k]],".rds")
  saveRDS(centered.delta.conf, file = file)
  

}

plan(sequential)


rm(list = ls())



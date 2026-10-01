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
source("weather temperature/trend analysis/functions.R")




bw_month = function(data.1, data.2){
  
  bw_list = matrix(NA, nrow = 3, ncol = 12,dimnames = list(c("period 1", "delta", "period 2"), month.name)) |> as.data.frame()
  
  results = future_sapply(c(1,4,7,10), function(m){
    
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    
    sample.2  = data.2 %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data.2)[2])
    sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
    K.2  = length(unique(sample.2$Year))
    
    sample.1 = data.1 %>% filter(MONTH %in% month.name[m.neighbor])|> dplyr::select(1,4:dim(data.1)[2])
    sample.1 = sample.1[rowSums(is.na(sample.1)) == 0,]
    K.1 = length(unique(sample.1$Year))
    
    bw.2      = k.fold.hc.cv(sample = sample.2, h.seq = seq(0.09, 0.16, 0.005),  deg = 2, K = K.2)  
    bw.delta  = k.fold.hc.cv(sample = sample.1, h.seq = seq(0.09, 0.16, 0.005),  deg = 2, diff = T, K = K.1, sample.2, bw.2) 
    bw.1      = k.fold.hc.cv(sample = sample.1, h.seq = seq(0.09, 0.16, 0.005),  deg = 2, K = K.1)  
    
    c(bw.2,bw.delta,bw.1)
    
  }, future.seed = T)
  
  bw_list[1, c(12,1:11)] = rep(results[1,],each = 3)
  bw_list[2, c(12,1:11)] = rep(results[2,],each = 3)
  bw_list[3, c(12,1:11)] = rep(results[3,],each = 3)
  
  return(bw_list)
}






##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
################                                            ##################################################
################   Trend analysis of time period 1952-1972  ##################################################
################                                            ##################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################




##############################################################################################################
###                                                    #######################################################
###   Mean and difference function: Cross Validation   #######################################################
###                                                    #######################################################
##############################################################################################################
                                                                                             #################
# Optional: For parallelizing computations                                                   #################
options(future.globals.maxSize = 64 * 1024^3)                                                #################
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                 #################
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                 #################
                                                                                             #################
for(k in c(1,2)){                                                                            #################
                                                                                             #################
  #                   1 (Berlin)                                                             #################
  data.example = list("Berlin", "Hamburg")                                                   #################
  station      = list("00433","01975")                                                       #################
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                  #################
                                                                                             #################  
                                                                                             #################
  if(k == 1){         s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1957 | Year == 1957  & MONTH   %in% month.name[1:4])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1957 | Year == 1957  & !(MONTH %in% month.name[1:4]))
    }else if(k == 2){ 
                      s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1968 | Year == 1968  & MONTH   %in% month.name[1:6])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1968 | Year == 1968  & !(MONTH %in% month.name[1:6]))
      }
  
                                                                                             #################
                                                                                             ################# 
  Bandwidths = bw_month(s.part1,s.part2)                                                          ############
  file = paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_bw_",data.example[[k]],".rds")
  saveRDS(Bandwidths, file)                                                                       ############
  print("save done")                                                                         #################
}                                                                                            #################
plan(sequential)                                                                             #################
                                                                                             #################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


##############################################################################################################
###                                                  #########################################################
###   Calculating long-run kernels (0,...,max.lag)   #########################################################
###                                                  #########################################################
##############################################################################################################
                                                                                                            ##
# Optional: For parallelizing computations                                                                  ##
options(future.globals.maxSize = 64 * 1024^3)                                                               ##
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                                ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                                ##
                                                                                                            ##
for(k in c(1,2)){                                                                                           ##   
#                   1 (Berlin)                                                                              ##
#                   2 (Hamburg)                                                                             ##
                                                                                                            ##
  data.example = list("Berlin","Hamburg")                                                                   ##
  station      = list("00433","01975")                                                                      ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                                 ##
                                                                                                            ##  
  # Seperate in two time periods                                                                            ##  
                                                                                                            ##
  if(k == 1){         s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1957 | Year == 1957  & MONTH   %in% month.name[1:4])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1957 | Year == 1957  & !(MONTH %in% month.name[1:4]))
    }else if(k == 2){ 
                      s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1968 | Year == 1968  & MONTH   %in% month.name[1:6])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1968 | Year == 1968  & !(MONTH %in% month.name[1:6]))
      }       
                                                                                                             
  # Maximum lag for each month                                                                              ##
  max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds"))##
                                                                                                            ##
  cov.1.list = list()                                                                                       ##
  cov.2.list = list()                                                                                       ##
                                                                                                            ##
                                                                                                            ##
  for(m in 1:12){                                                                                           ## 
                                                                                                            ##
  part1  = s.part1 |>                                                                                       ##
            filter(MONTH == month.name[m]) |>                                                               ##
              dplyr::select(1,4:dim(s.part1)[2])                                                            ##
                                                                                                            ##
  p1 = dim(part1)[2] - 1                                                                                    ##
                                                                                                            ##
  part2  = s.part2 |>                                                                                       ##
            filter(MONTH == month.name[m]) |>                                                               ##
              dplyr::select(1,4:dim(s.part2)[2])                                                            ##
                                                                                                            ##
  p2 = dim(part2)[2] - 1                                                                                    ##
                                                                                                            ##
  lag.k.Gamma.1 = list()                                                                                    ##
  lag.k.Gamma.2 = list()                                                                                    ##
  
  file = paste0("weather temperature/kernel weights/full_w_s_lag0.rds") 
  w   = readRDS(file)                                                                              
  
  lag.Gamma.1   = eval.weights(w,observation.transformation(part1[,-1],grid.type = "less", periodic = T, m = 24))
  lag.k.Gamma.1[[1]] = lag.Gamma.1
  
  lag.Gamma.2   = eval.weights(w,observation.transformation(part2[,-1],grid.type = "less", periodic = T, m = 24))
  lag.k.Gamma.2[[1]] = lag.Gamma.2
  
  rm(w, lag.Gamma.2,lag.Gamma.1)
  
  lag.k.Gamma.part1 = future_lapply(1:max.lag[m,2], function(k) {
    
    file  = paste0("weather temperature/kernel weights/full_w_s_lag",k,".rds")
    w.lag = readRDS(file)
    n.year = unique(part1$Year) 
    lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(part1[which(part1$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 24), lag = k)}))/length(n.year)
    lag.Gamma
  })
  lag.k.Gamma.1 = c(lag.k.Gamma.1, lag.k.Gamma.part1)
  
  lag.k.Gamma.part2 = future_lapply(1:max.lag[m,2], function(k) { 
    
    file  = paste0("weather temperature/kernel weights/full_w_s_lag",k,".rds")
    w.lag = readRDS(file)
    n.year = unique(part2$Year) 
    lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(part2[which(part2$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 24), lag = k)}))/length(n.year)
    lag.Gamma
  })                                                                                                        ##
  lag.k.Gamma.2 = c(lag.k.Gamma.2,lag.k.Gamma.part2)                                                        ##
                                                                                                            ##
  cov.1.list[[m]]   = lag.k.Gamma.1                                                                         ##
  cov.2.list[[m]]   = lag.k.Gamma.2                                                                         ## 
                                                                                                            ##
  rm(lag.k.Gamma.1,p1)                                                                                      ##
  rm(lag.k.Gamma.2,p2)                                                                                      ##
                                                                                                            ##
  rm(part1, part2, lag.k.Gamma.s.part1, lag.k.Gamma.part2)                                                  ##
  print(paste0("done: ", month.name[m]))                                                                    ##
  }                                                                                                         ##
                              

  saveRDS(cov.1.list , paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_list_Gamma_period1_",data.example[[k]],".rds"))
  saveRDS(cov.2.list , paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_list_Gamma_period2_",data.example[[k]],".rds"))
    

                                                                                                            ##
}                                                                                                           ##
plan(sequential)                                                                                            ##
                                                                                                            ##
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

















###############################################################################################################
###                                                                                     #######################
###   Constructing confidence intervals and bands (integral and non centered process)   #######################
###                                                                                     #######################
###############################################################################################################
                                                                                                             ##
# Optional: For parallelizing computations                                                                   ##
options(future.globals.maxSize = 64 * 1024^3)                                                                ##
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                                 ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                                 ##
                                                                                                             ##
for(k in c(1,2)){                                                                                            ##   
#                   1 (Berlin)                                                                               ##
#                   2 (Hamburg)                                                                              ##
                                                                                                             ##   
  data.example     = list("Berlin","Hamburg")                                                                ##
  station          = list("00433","01975")                                                                   ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                                  ##
                                                                                                             ##
  if(k == 1){         s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1957 | Year == 1957  & MONTH   %in% month.name[1:4])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1957 | Year == 1957  & !(MONTH %in% month.name[1:4]))
    }else if(k == 2){ 
                      s.part1 = data.s.34h |>                                                                 
                                filter(Year < 1968 | Year == 1968  & MONTH   %in% month.name[1:6])                                            
                      s.part2 = data.s.34h |>                                                                 
                                filter(Year > 1968 | Year == 1968  & !(MONTH %in% month.name[1:6]))
      } 
      
  
  # Mean and difference function: Estimation                                                                                              
  Bandwidths = readRDS(paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_bw_",data.example[[k]],".rds"))
  est        = est.results(s.part1, s.part2, Bandwidths,from = 1,to = 12)                                   
                                                                                                         

  # Loading estimated lag covariance kernels for lag = 0,...,12 
  cov.1.list = readRDS(paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_list_Gamma_period1_",data.example[[k]],".rds")) 
  cov.2.list = readRDS(paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_list_Gamma_period2_",data.example[[k]],".rds")) 
  
  # Loading maximum lag of long run kernel 
  max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds")) ## 
                                                                                                             ##
                                                                                                             ## 
  start.24 = 21;end.24 = 117                                                                                 ##
                                                                                                             ##
############################################################################################################### 
## Calculating variance for each month ########################################################################
###############################################################################################################        
                                                                                                             ##
  lr.var.1.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,                                         ##
                                                                lapply(1:(max.lag[m,2]+1),  
                                                                       function(j){if(j >= 2){
                                                                         (cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                       }else{cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
  }))
  
  
  lr.var.2.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,
                                                              lapply(1:(max.lag[m,2]+1),  
                                                                     function(j){if(j >= 2){
                                                                       (cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                     }else{cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
  }))
                                                                                                             ##
                                                                                                             ##
###############################################################################################################
### Integral of delta test with construction of confidence bands ##############################################
###############################################################################################################                  
                                                                                                             ##
  set.seed(2005)
  args.MB       = list(est = est, bandwidth = Bandwidths, B = 10000, dependent = T)
  integral.test = test.int(s.part1, s.part2, unique(est$delta_int$ESTIMATE), lr.var.1.month, lr.var.2.month, from = 1, to = 12, test = "two-sided", alpha = rep(0.95,12), approx = "MB", args.MB)
  
  file = paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_quantile_integral_",data.example[[k]],".rds")
  saveRDS(integral.test, file = file)

  integral.conf       = tibble(TIME = rep(hms::as_hms(c(as.POSIXct("1970-01-01 00:00:00"),as.POSIXct("1970-01-01 23:59:59"))),times = 12))
  integral.conf$UP    = rep(unlist(integral.test$confInterval)[seq(2,24,2)], each = 2)
  integral.conf$LO    = rep(unlist(integral.test$confInterval)[seq(1,24,2)], each = 2)
  integral.conf$MONTH =  factor(rep(month.name , each = 2), level = month.name)

  file = paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_CB_integral_",data.example[[k]],".rds")#
  saveRDS(integral.conf, file = file)                                                                        
                                                                                                             ## 
                                                                                                             ##
                                                                                                             ##
                                                                                                             ##
                                                                                                             ##  
###############################################################################################################            
### Calculating diagonal of long run covariance matrix for each month #########################################
###############################################################################################################            
                                                                                                             ##
  lr.Gamma.1.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                           function(j){if(j >= 2){
                                                                             (cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                           }else{cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                           }))
  return(diag(lr.cov))})
  
  lr.Gamma.2.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                           function(j){if(j >= 2){
                                                                             (cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                           }else{cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                           }))
  return(diag(lr.cov))})
                                                                                                            

  lr.Gamma = lapply(1:12, function(m) { sample.2  = s.part2 %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(s.part2)[2])
                                        sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
  
                                        sample.1  = s.part1 %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(s.part1)[2])
                                        sample.1  = sample.1[rowSums(is.na(sample.1)) == 0,]
  
                                        return(lr.Gamma.1.month[[m]] + dim(sample.1)[1] / dim(sample.2)[1] * lr.Gamma.2.month[[m]])
  }) 

    
                                                                                                            ##
##############################################################################################################               
### Dependent Multiplier Bootstrap on non centered difference function for each month ########################   
##############################################################################################################            
                                                                                                            ## 
  set.seed(2005)                                                                                            ##
  q.list = q.month(s.part1, s.part2, Bandwidths, est, lr.Gamma,                                             ##
                   from = 1, to = 12, alpha = rep(0.95,12), B = 10000, depend = T, int = F, grid = "sparse")##
                                                                                                            ##
  file = paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_quantile_",data.example[[k]],".rds") ##
  saveRDS(q.list, file = file)                                                                              ##          
                                                                                                            ##
  # Skip:                                                                                                   
  #q.list = readRDS(paste0("weather temperature/trend analysis/Supplement Results part2/sparse_quantile_",data.example[[k]],".rds"))
  
  delta.conf = CB(s.part1, est, lr.Gamma, q.list[1,])                                                       ##
  file = paste0("weather temperature/trend analysis/relocation results/sparse_",station[[k]],"_CB_delta_",data.example[[k]],".rds") ##
  saveRDS(delta.conf, file = file)                                                                          ##
                                                                                                            ##
}                                                                                                           ##
                                                                                                            ## 
plan(sequential)                                                                                            ##                                                                             
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


  
  
  





















bw_month = function(data.1, data.2){
  
  bw_list = matrix(NA, nrow = 3, ncol = 12,dimnames = list(c("period (2014-11 - 2025)", "delta", "period (2000 - 2014-10)"), month.name)) |> as.data.frame()
  
  results = future_sapply(c(1,4,7,10), function(m){
    
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    
    sample.2  = data.2 %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data.2)[2])
    sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
    K.2  = length(unique(sample.2$Year))
    
    sample.1 = data.1 %>% filter(MONTH %in% month.name[m.neighbor])|> dplyr::select(1,4:dim(data.1)[2])
    sample.1 = sample.1[rowSums(is.na(sample.1)) == 0,]
    K.1 = length(unique(sample.1$Year))
    
    bw.2      = k.fold.hc.cv(sample = sample.2, h.seq = seq(0.045, 0.1,  0.001),  deg = 2, K = K.2)  
    bw.delta  = k.fold.hc.cv(sample = sample.1, h.seq = seq(0.055, 0.15, 0.001),  deg = 2, diff = T, K = K.1, sample.2, bw.2) 
    bw.1      = k.fold.hc.cv(sample = sample.1, h.seq = seq(0.045, 0.1,  0.001),  deg = 2, K = K.1)  
    
    c(bw.2,bw.delta,bw.1)
    
  }, future.seed = T)
  
  bw_list[1, c(12,1:11)] = rep(results[1,],each = 3)
  bw_list[2, c(12,1:11)] = rep(results[2,],each = 3)
  bw_list[3, c(12,1:11)] = rep(results[3,],each = 3)
  
  return(bw_list)
}






##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
################                                            ##################################################
################   Trend analysis of time period 2000-2025  ##################################################
################                                            ##################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################




##############################################################################################################
###                                                    #######################################################
###   Mean and difference function: Cross Validation   #######################################################
###                                                    #######################################################
##############################################################################################################
                                                                                             #################
# Optional: For parallelizing computations                                                   #################
options(future.globals.maxSize = 64 * 1024^3)                                                #################
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                 #################
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                 #################
                                                                                             #################
                                                                                             #################
  # 1 (Frankfurt am Main)                                                                    #################
  k = 1                                                                                      #################
  data.example = list( "Frankfurt_Main" )                                                    #################
  station      = list("01420")                                                               #################
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                  #################                                                                                           #################  
  # Seperate in two time periods                                                             #################  
                                                                                             #################
                                                                                             #################
  d.part1 = data.d.34h |> filter(Year < 2014 | Year == 2014  &   MONTH %in% month.name[1:10])   ##############                           
  d.part2 = data.d.34h |> filter(Year > 2014 | Year == 2014  & !(MONTH %in% month.name[1:10]))  ##############
                                                                                             #################
                                                                                             #################
                                                                                             ################# 
  Bandwidths = bw_month(d.part1,d.part2)                                                         #############
  file = paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_bw_",data.example[[k]],".rds") 
  saveRDS(Bandwidths, file)                                                                      #############
  print("save done")                                                                         #################
                                                                                             #################
plan(sequential)                                                                             #################
                                                                                             #################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################



##############################################################################################################
###                                                  #########################################################
###   Calculating long-run kernels (0,...,max.lag)   #########################################################
###                                                  #########################################################
##############################################################################################################
                                                                                                            ##
# Optional: For parallelizing computations                                                                  ##
options(future.globals.maxSize = 64 * 1024^3)                                                               ##
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                                ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                                ##
                                                                                                            ##
                                                                                                            ##
  # 1 (Frankfurt am Main)                                                                                   ##
  k = 1                                                                                                     ##
  data.example = list( "Frankfurt_Main" )                                                                   ##
  station      = list("01420")                                                                              ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                                 ##
                                                                                                            ##
  # Seperate in two time periods: (2000 - 2014-10) and (2014-11 - 2025)                                     ##           ##  
                                                                                                            ##
  d.part1 = data.d.34h |> filter(Year < 2014 | Year == 2014  &   MONTH %in% month.name[1:10])               ##                             
  d.part2 = data.d.34h |> filter(Year > 2014 | Year == 2014  & !(MONTH %in% month.name[1:10]))              ##
                                                                                                            ##
  # Maximum lag for each month                                                                              ##
  max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds"))##
                                                                                                            ##
  cov.1.list = list()                                                                                       ##
  cov.2.list = list()                                                                                       ##
                                                                                                            ##
  file = paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[k]],".rds")     ##
  bandwidth.d = readRDS(file)                                                                               ##

  for(m in 1:12){                                                                                             ## 
                                                                                                              ##
    part1  = d.part1 |>                                                                                       ##
              filter(MONTH == month.name[m]) |>                                                               ##
                dplyr::select(1,4:dim(d.part1)[2])                                                            ##
                                                                                                              ##
    p1 = dim(part1)[2] - 1                                                                                    ##
                                                                                                              ##
    part2  = d.part2 |>                                                                                       ##
              filter(MONTH == month.name[m]) |>                                                               ##
                dplyr::select(1,4:dim(d.part2)[2])                                                            ##
                                                                                                              ##
    p2 = dim(part2)[2] - 1                                                                                    ##
                                                                                                              ##
    lag.k.Gamma.1 = list()                                                                                    ##
    lag.k.Gamma.2 = list()                                                                                    ##
    
    if(bandwidth.d[m] < 0.1){ file.d = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%03d", as.integer(bandwidth.d[m] * 100)),".rds")
    }else{                    file.d = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%02d", as.integer(bandwidth.d[m] * 10)),".rds")}
    
    wd  = readRDS(file.d)                                                                               
    
    lag.Gamma.1   = eval.weights(wd,observation.transformation(part1[,-1],grid.type = "less", periodic = T, m = 144))
    lag.k.Gamma.1[[1]] = lag.Gamma.1
    
    lag.Gamma.2   = eval.weights(wd,observation.transformation(part2[,-1],grid.type = "less", periodic = T, m = 144))
    lag.k.Gamma.2[[1]] = lag.Gamma.2
    
    rm(wd, lag.Gamma.2,lag.Gamma.1)
    
    lag.k.Gamma.part1 = future_lapply(1:max.lag[m,2], function(k) {
      
      if(bandwidth.d[m] < 0.1){ file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%03d", as.integer(bandwidth.d[m] * 100)),".rds")}
      else{                     file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%02d", as.integer(bandwidth.d[m] * 10)),".rds")}
      w.lag = readRDS(file.d)
      n.year = unique(part1$Year) 
      lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(part1[which(part1$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 144), lag = k)}))/length(n.year)
      lag.Gamma
    })
    lag.k.Gamma.1 = c(lag.k.Gamma.1, lag.k.Gamma.part1)
    
    lag.k.Gamma.part2 = future_lapply(1:max.lag[m,2], function(k) { 
      
      if(bandwidth.d[m] < 0.1){ file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%03d", as.integer(bandwidth.d[m] * 100)),".rds")}
      else{                     file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%02d", as.integer(bandwidth.d[m] * 10)),".rds")}
      w.lag = readRDS(file.d)
      n.year = unique(part2$Year) 
      lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(part2[which(part2$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 144), lag = k)}))/length(n.year)
      lag.Gamma
    })                                                                                                      ##
    lag.k.Gamma.2 = c(lag.k.Gamma.2,lag.k.Gamma.part2)                                                      ##
                                                                                                            ##
    cov.1.list[[m]]   = lag.k.Gamma.1                                                                       ##
    cov.2.list[[m]]   = lag.k.Gamma.2                                                                       ## 
                                                                                                            ##
    rm(lag.k.Gamma.1,p1)                                                                                    ##
    rm(lag.k.Gamma.2,p2)                                                                                    ##
                                                                                                            ##
    rm(part1, part2, lag.k.Gamma.d.part1, lag.k.Gamma.part2)                                                ##
    print(paste0("done: ", month.name[m]))                                                                  ##
  }                                                                                                         ##
  
  saveRDS(cov.1.list , paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_list_Gamma_period1_",data.example[[k]],".rds"))
  saveRDS(cov.2.list , paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_list_Gamma_period2_",data.example[[k]],".rds"))
  
                                                                                                            ##
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

















###############################################################################################################
###                                                                                     #######################
###   Constructing confidence intervals and bands (integral and non centered process)   #######################
###                                                                                     #######################
###############################################################################################################
                                                                                                             ##
# Optional: For parallelizing computations                                                                   ##
options(future.globals.maxSize = 64 * 1024^3)                                                                ##
plan(multisession, workers = 120) # MaRC3a: partition=mqtest                                                 ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal                                                 ##
                                                                                                             ##
  # 1 (Frankfurt am Main)                                                                                    ##
  k = 1                                                                                                      ##                                                                                                          ##   
  data.example = list( "Frankfurt_Main" )                                                                    ##
  station      = list("01420")                                                                               ##
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))                                  ##
                                                                                                             ##
  d.part1 = data.d.34h |> filter(Year < 2014 | Year == 2014  &   MONTH %in% month.name[1:10])                                            
  d.part2 = data.d.34h |> filter(Year > 2014 | Year == 2014  & !(MONTH %in% month.name[1:10]))
  
  
  # Mean and difference function: Estimation                                                                                              
  Bandwidths = readRDS(paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_bw_",data.example[[k]],".rds"))   
  est        = est.results(d.part1, d.part2, Bandwidths,from = 1,to = 12)                                   
  
  
  # Loading estimated lag covariance kernels for lag = 0,...,12 
  cov.1.list = readRDS(paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_list_Gamma_period1_",data.example[[k]],".rds")) 
  cov.2.list = readRDS(paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_list_Gamma_period2_",data.example[[k]],".rds")) 
  
  # Loading maximum lag of long run kernel 
  max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds"))  ## 
                                                                                                              ##
                                                                                                              ## 
  start.24 = 21;end.24 = 117                                                                                  ##
                                                                                                              ##
  ##############################################################################################################
  ## Calculating variance for each month #######################################################################
  ##############################################################################################################        
                                                                                                              ##
  lr.var.1.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,                                          ##
                                                                lapply(1:(max.lag[m,2]+1),  
                                                                       function(j){if(j >= 2){
                                                                         (cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                       }else{cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
  }))
  
  
  lr.var.2.month = unlist(lapply(1:12, function(m) {mean(Reduce(`+`,
                                                                lapply(1:(max.lag[m,2]+1),  
                                                                       function(j){if(j >= 2){
                                                                         (cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                       }else{cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]} })))
  }))
                                                                                                               ##
                                                                                                               ##
  ###############################################################################################################
  ### Integral of delta test with construction of confidence bands ##############################################
  ###############################################################################################################                  
                                                                                                               ##
  set.seed(2005)
  args.MB       = list(est = est, bandwidth = Bandwidths, B = 10000, dependent = T)
  integral.test = test.int(d.part1, d.part2, unique(est$delta_int$ESTIMATE), lr.var.1.month, lr.var.2.month, from = 1, to = 12, test = "two-sided", alpha = rep(0.95,12), approx = "MB", args.MB)
  
  file = paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_quantile_integral_",data.example[[k]],".rds")
  saveRDS(integral.test, file = file)
  
  integral.conf       = tibble(TIME = rep(hms::as_hms(c(as.POSIXct("1970-01-01 00:00:00"),as.POSIXct("1970-01-01 23:59:59"))),times = 12))
  integral.conf$UP    = rep(unlist(integral.test$confInterval)[seq(2,24,2)], each = 2)
  integral.conf$LO    = rep(unlist(integral.test$confInterval)[seq(1,24,2)], each = 2)
  integral.conf$MONTH =  factor(rep(month.name , each = 2), level = month.name)
  
  file = paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_CB_integral_",data.example[[k]],".rds")     
  saveRDS(integral.conf, file = file)                                                                          ##
                                                                                                               ## 
                                                                                                               ##
                                                                                                               ##
                                                                                                               ##
                                                                                                               ##  
  ###############################################################################################################            
  ### Calculating diagonal of long run covariance matrix for each month #########################################
  ###############################################################################################################            
                                                                                                               ##
  lr.Gamma.1.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                           function(j){if(j >= 2){
                                                                             (cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                           }else{cov.1.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                           }))
  return(diag(lr.cov))})
  
  lr.Gamma.2.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                           function(j){if(j >= 2){
                                                                             (cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                           }else{cov.2.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                           }))
  return(diag(lr.cov))})
  
  lr.Gamma = lapply(1:12, function(m) { sample.2  = d.part2 %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(d.part2)[2])
                                        sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
  
                                        sample.1  = d.part1 %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(d.part1)[2])
                                        sample.1  = sample.1[rowSums(is.na(sample.1)) == 0,]
  
                                        return(lr.Gamma.1.month[[m]] + dim(sample.1)[1] / dim(sample.2)[1] * lr.Gamma.2.month[[m]])
  })
  
  ############################################################################################################             
  ### Dependent Multiplier Bootstrap on non centered difference function for each month ######################
  ############################################################################################################          
                                                                                                            ## 
  set.seed(2005)                                                                                            ##
  q.list = q.month(d.part1, d.part2, Bandwidths, est, lr.Gamma,                                             ##
                   from = 1, to = 12, alpha = rep(0.95,12), B = 10000, depend = T, int = F, grid = "dense") ##
                                                                                                            ##
  file = paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_quantile_",data.example[[k]],".rds")     
  saveRDS(q.list, file = file)                                                                              ##
                                                                                                            ##
  # Skip:                                                                                                   
  #q.list = readRDS(paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_quantile_",data.example[[k]],".rds"))                                                                                                      
  
  delta.conf = CB(d.part1, est, lr.Gamma, q.list[1,])                                                       ##
  file = paste0("weather temperature/trend analysis/relocation results/dense_",station[[k]],"_CB_delta_",data.example[[k]],".rds")      
  saveRDS(delta.conf, file = file)                                                                          ##
                                                                                                            ##
                                                                                                           ##
plan(sequential)                                                                                            ##    
                                                                                                            ##
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
  
  
  
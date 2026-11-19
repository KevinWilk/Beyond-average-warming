library(ggplot2)
library(reshape2)
library(locpol)
library(biLocPol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(tidyr)
library(lubridate)
library(hms)
library(dplyr)
library(tibble)
library(plotly)
library(gridExtra)
library(CompQuadForm)
library(tensorA)



source("Mod_biLocPol.R")
source("weather temperature/functions.R")
source("weather temperature/long run kernel/function_test.R")

############################################################################
k = 1  # Change to: 1 (Berlin)                                            ##
#                   2 (Frankfurt am Main)                                 ##
#                   3 (Hamburg)                                           ##
#                   4 (Munich)                                            ##
data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      ##
load(paste0("weather temperature/data sets/",data.example[[k]],".RData")) ##
############################################################################



########################################################
# Estimating lag covariance kernels for lag = 0,...,12 #
########################################################


cov.s.list = list()
cov.d.list = list()


# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 64 * 1024^3)                ##   
plan(multisession, workers = 120) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################

bandwidth.d = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_Gamma_d_",data.example[[k]],".rds"))

for(m in 1:12){
  
    sample.dense  = data.d.34h |> 
                      filter(MONTH == month.name[m]) |> 
                        dplyr::select(1,4:dim(data.d.34h)[2])
    
    sample.dense = sample.dense[rowSums(is.na(sample.dense)) == 0,] 
    pd = dim(sample.dense)[2] - 1
  
    sample.sparse = data.s.34h |> 
                      filter(MONTH == month.name[m]) |> 
                        dplyr::select(1,4:dim(data.s.34h)[2])
    
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]  
    p = dim(sample.sparse)[2] - 1
  
    lag.k.Gamma   = list()
    lag.k.Gamma.d = list()  
  
    file.s = paste0("weather temperature/kernel weights/full_w_s_lag0.rds")
    if(bandwidth.d[m] < 0.1){ file.d = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%03d", as.integer(bandwidth.d[m] * 100)),".rds")
    }else{                    file.d = paste0("weather temperature/kernel weights/full_w_d_lag0_",sprintf("%02d", as.integer(bandwidth.d[m] * 10)),".rds")}
    w   = readRDS(file.s)
    wd  = readRDS(file.d)
  
    lag.Gamma = eval.weights(w,observation.transformation(sample.sparse[,-1],grid.type = "less", periodic = T, m = 24))
    lag.k.Gamma[[1]] = lag.Gamma

    lag.Gamma.d   = eval.weights(wd,observation.transformation(sample.dense[,-1],grid.type = "less", periodic = T, m = 144))
    lag.k.Gamma.d[[1]] = lag.Gamma.d
  
    rm(w,lag.Gamma)
    rm(wd,lag.Gamma.d)
  
    lag.k.Gamma.part2 = future_lapply(1:12, function(k) {
      file.s = paste0("weather temperature/kernel weights/full_w_s_lag",k,".rds")
      w.lag = readRDS(file.s)
      n.year = unique(sample.sparse$Year) 
      lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(sample.sparse[which(sample.sparse$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 24), lag = k)}))/length(n.year)
      lag.Gamma
    })
    lag.k.Gamma = c(lag.k.Gamma,lag.k.Gamma.part2)
  
    lag.k.Gamma.d.part2 = future_lapply(1:12, function(k) {  
      if(bandwidth.d[m] < 0.1){ file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%03d", as.integer(bandwidth.d[m] * 100)),".rds")}
      else{                     file.d = paste0("weather temperature/kernel weights/full_w_d_lag",k,"_",sprintf("%02d", as.integer(bandwidth.d[m] * 10)),".rds")}
      w.lag = readRDS(file.d)
      n.year = unique(sample.dense$Year) 
      lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(sample.dense[which(sample.dense$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 144), lag = k)}))/length(n.year)
      lag.Gamma
    })
    lag.k.Gamma.d = c(lag.k.Gamma.d,lag.k.Gamma.d.part2)
  
    cov.s.list[[m]]   = lag.k.Gamma
    cov.d.list[[m]]   = lag.k.Gamma.d
  
    rm(lag.k.Gamma,p)
    rm(lag.k.Gamma.d,pd)
  
    rm(sample.dense,sample.sparse,lag.k.Gamma.d.part2,lag.k.Gamma.part2)
    print(paste0("done: ", month.name[m]))
}


saveRDS(cov.d.list , paste0("weather temperature/long run kernel/Results/list_Gamma_d_",data.example[[k]],".rds"))
saveRDS(cov.s.list , paste0("weather temperature/long run kernel/Results/list_Gamma_s_",data.example[[k]],".rds"))



#############################################################################################################################
# Skip: Estimating lag covariance kernels for lag = 0,...,12 ################################################################
# cov.s.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_s_",data.example[[k]],".rds")) ########
# cov.d.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_d_",data.example[[k]],".rds")) ########
#############################################################################################################################







#######################
# Finding maximal lag #
#######################


start.24 = 21; end.24 = 117 

hat.rho = function(h, cov.data){lapply(1:12, function(m) {cov.lag    = sqrt(mean((cov.data[[m]][[h]][start.24:end.24,start.24:end.24])^2))
                                                          cov        = mean(diag( cov.data[[m]][[1]][start.24:end.24,start.24:end.24]))
                                                          cov.lag/cov })}

rho.hat.s.df = as.data.frame(do.call(rbind, lapply(2:13, function(k) {unlist(hat.rho(k, cov.s.list))})))
colnames(rho.hat.s.df) = month.name
rownames(rho.hat.s.df) = paste0("lag ", 1:12)
rho.hat.s.df = rho.hat.s.df |> rownames_to_column("lag") |> mutate(lag = as.integer(gsub("lag\\s*", "", lag))) |> pivot_longer(cols = -lag,names_to = "MONTH",values_to = "value") |> mutate(MONTH = factor(MONTH, levels = month.name))

rho.hat.d.df = as.data.frame(do.call(rbind, lapply(2:13, function(k) {unlist(hat.rho(k, cov.d.list))})))
colnames(rho.hat.d.df) = month.name
rownames(rho.hat.d.df) = paste0("lag ", 1:12)
rho.hat.d.df = rho.hat.d.df |> rownames_to_column("lag") |> mutate(lag = as.integer(gsub("lag\\s*", "", lag))) |> pivot_longer(cols = -lag,names_to = "MONTH",values_to = "value") |> mutate(MONTH = factor(MONTH, levels = month.name))


test.s.val = lapply(1:12, function(month) {inference.autocovariance.test(data.s.24h   |> select(1:27) , m = month, alpha = 0.95, max.lag = 12)})
test.d.val = lapply(1:12, function(month) {inference.autocovariance.test(data.d.24h.hourly |> select(1:27) , m = month, alpha = 0.95, max.lag = 12)}) 

test.s.val = data.frame(lag = rep(c(0.6,2:11,12.5), times = 12), value = unlist(test.s.val), MONTH = factor(rep(month.name[1:length(test.s.val)], lengths(test.s.val)), levels = month.name))
test.d.val = data.frame(lag = rep(c(0.6,2:11,12.5), times = 12), value = unlist(test.d.val), MONTH = factor(rep(month.name[1:length(test.d.val)], lengths(test.d.val)), levels = month.name))


saveRDS(test.s.val, paste0("weather temperature/long run kernel/Results/test_max_lag_s_",data.example[[k]],".rds"))
saveRDS(test.d.val, paste0("weather temperature/long run kernel/Results/test_max_lag_d_",data.example[[k]],".rds"))


plan(sequential)


#############################################################################################################################
# Skip running test.s.val and test.d.val                     ################################################################
# test.s.val = readRDS(paste0("weather temperature/long run kernel/Results/test_max_lag_s_",data.example[[k]],".rds")) ######
# test.d.val = readRDS(paste0("weather temperature/long run kernel/Results/test_max_lag_d_",data.example[[k]],".rds")) ######
#############################################################################################################################




max.lag.d = rho.hat.d.df |>
              arrange(MONTH, lag) |>
                mutate(threshold = test.d.val$value) |>
                  group_by(MONTH) |>
                    summarise(first_lag = lag[which.max(value < threshold)[1]] - 1,.groups = "drop")

max.lag.s = rho.hat.s.df |>
              arrange(MONTH, lag) |>
                mutate(threshold = test.s.val$value) |>
                  group_by(MONTH) |>
                    summarise(first_lag = lag[which(value < threshold)[1]] - 1,.groups = "drop") 



max.lag = data.frame(max.lag.d$MONTH, max.lag = apply(cbind(max.lag.d$first_lag,max.lag.s$first_lag),1,max))

saveRDS(max.lag , paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds"))







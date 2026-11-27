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


##################################
### Global test over 12 months ###
##################################


for(k in 1:4){
  

  data.example     = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")     
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))



  #########################################################################################################################################
  #(1952-1972) with (2000-2025) ###########################################################################################################

  #Integral: approximation with normaldistribution
  #test.int = readRDS(paste0("weather temperature/confidence bands/Results/quantile_integral_",data.example[[k]],".rds"))$pval

  #Integral: approximation with MB
  test.int    = readRDS(paste0("weather temperature/confidence bands/Results/quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.int    = p.adjust(p.val, method = "holm")
  test.int    = format(test.int, scientific  = F)
  cat(paste0("\n Global integral difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))
  print(test.int)
  
  
  
  # not centered difference process
  test.delta  = readRDS(paste0("weather temperature/confidence bands/Results/quantile_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = est. integral"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.delta  = p.adjust(p.val, method = "holm")
  test.delta  = format(test.delta, scientific  = F)
  cat(paste0("\n Global difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))
  print(test.delta)
  
  
  
  # centered difference process
  test.center.delta = readRDS(paste0("weather temperature/confidence bands/Results/quantile_P_",data.example[[k]],".rds"))
  p.val             = data.frame(t(test.center.delta[2,]), row.names = c("H0: centered delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000                                                                                 
  test.center.delta = p.adjust(p.val, method = "holm")
  test.center.delta = format(test.center.delta, scientific  = F)
  cat(paste0("\n Global centered difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))
  print(test.center.delta)
 
  
  
  
  # Summary of all test: (1952-1972) with (2000-2025)
  years_2000_2025_with_1952_1972 = rbind(test.int,test.delta,test.center.delta)
  
  
  
  
  
  #########################################################################################################################################
  #                                                                  ######################################################################
  # Robustness check: Multiply bandwidth of delta estimator by 1.5   ######################################################################
  #                                                                  ######################################################################
  #Integral: approximation with MB ########################################################################################################                 
  test.int    = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 1.5/quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)                                  ###
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.int    = p.adjust(p.val, method = "holm")                                                                                        ###
  test.int    = format(test.int, scientific  = F)                                                                                       ###
  cat(paste0("\n Global integral difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))             ###
  print(test.int)                                                                                                                       ###
                                                                                                                                        ###
  # not centered difference process #######################################################################################################
  test.delta  = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 1.5/quantile_",data.example[[k]],".rds"))      ###
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = est. integral"))|> setNames(month.name)                        ###
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.delta  = p.adjust(p.val, method = "holm")                                                                                        ###
  test.delta  = format(test.delta, scientific  = F)                                                                                     ###
  cat(paste0("\n Global difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))                      ###
  print(test.delta)                                                                                                                     ###
                                                                                                                                        ###
  # centered difference process                                                                                                         ###
  test.center.delta = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 1.5/quantile_P_",data.example[[k]],".rds"))#
  p.val             = data.frame(t(test.center.delta[2,]), row.names = c("H0: centered delta = 0"))|> setNames(month.name)              ###
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.center.delta = p.adjust(p.val, method = "holm")                                                                                  ###
  test.center.delta = format(test.center.delta, scientific  = F)                                                                        ### 
  cat(paste0("\n Global centered difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))             ###
  print(test.center.delta)                                                                                                              ###
  #########################################################################################################################################
  #########################################################################################################################################
  
  robustness_15_check_years_2000_2025_with_1952_1972 = rbind(test.int,test.delta,test.center.delta)
  
  
  #########################################################################################################################################
  #                                                                  ######################################################################
  # Robustness check: Multiply bandwidth of delta estimator by 0.75  ###################################################################### 
  #                                                                  ######################################################################
  #Integral: approximation with MB ########################################################################################################                 
  test.int    = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 0.75/quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)                                  ###
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.int    = p.adjust(p.val, method = "holm")                                                                                        ###
  test.int    = format(test.int, scientific  = F)                                                                                       ###
  cat(paste0("\n Global integral difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))             ###
  print(test.int)                                                                                                                       ###
                                                                                                                                        ###
  # not centered difference process #######################################################################################################
  test.delta  = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 0.75/quantile_",data.example[[k]],".rds"))     ###
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = est. integral"))|> setNames(month.name)                        ###
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.delta  = p.adjust(p.val, method = "holm")                                                                                        ###
  test.delta  = format(test.delta, scientific  = F)                                                                                     ###
  cat(paste0("\n Global difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))                      ###
  print(test.delta)                                                                                                                     ###
                                                                                                                                        ###
  # centered difference process                                                                                                         ###
  test.center.delta = readRDS(paste0("weather temperature/confidence bands/Results/multiply by 0.75/quantile_P_",data.example[[k]],".rds"))
  p.val             = data.frame(t(test.center.delta[2,]), row.names = c("H0: centered delta = 0"))|> setNames(month.name)              ###             
  p.val[p.val == 0] = 1/10000                                                                                                           ###
  test.center.delta = p.adjust(p.val, method = "holm")                                                                                  ###
  test.center.delta = format(test.center.delta, scientific  = F)                                                                        ### 
  cat(paste0("\n Global centered difference test of (1952-1972) and (2000-2025) with multiplier bootstrap approximation:"))             ###
  print(test.center.delta)                                                                                                              ###
  #########################################################################################################################################
  #########################################################################################################################################

  robustness_075_check_years_2000_2025_with_1952_1972 = rbind(test.int,test.delta,test.center.delta)
  
  
  
  
  
  
  
  
  
  #########################################################################################################################################
  #(1952-1962) with (1963-1972) ###########################################################################################################
  
  #Integral: approximation with MB
  test.int = readRDS(paste0("weather temperature/trend analysis/trend results/sparse_quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.int    = p.adjust(p.val, method = "holm")
  test.int    = format(test.int, scientific  = F)
  cat(paste0("\n Global integral difference test of (1952-1962) and (1963-1972) with multiplier bootstrap approximation:"))
  print(test.int)
  
  
  
  # not centered difference process
  test.delta  = readRDS(paste0("weather temperature/trend analysis/trend results/sparse_quantile_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.delta  = p.adjust(p.val, method = "holm")
  test.delta  = format(test.delta, scientific  = F)
  cat(paste0("\n Global difference test of (1952-1962) and (1963-1972) with multiplier bootstrap approximation:"))
  print(test.delta)
  
  
  
  # centered difference process
  test.center.delta = readRDS(paste0("weather temperature/trend analysis/trend results/sparse_P_quantile_",data.example[[k]],".rds"))
  p.val             = data.frame(t(test.center.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.center.delta = p.adjust(p.val, method = "holm")
  test.center.delta = format(test.center.delta, scientific  = F)
  cat(paste0("\n Global difference test of (1952-1962) and (1963-1972) with multiplier bootstrap approximation:"))
  print(test.center.delta)
  
  
  
  # Summary of all test: (1952-1962) with (1963-1972)
  years_1952_1962_with_1963_1972 = rbind(test.int,test.delta,test.center.delta)
  
  #########################################################################################################################################
  #(2000-2012) with (2013-2025) ###########################################################################################################
  
  #Integral: approximation with MB
  test.int    = readRDS(paste0("weather temperature/trend analysis/trend results/dense_quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.int    = p.adjust(p.val, method = "holm")
  test.int    = format(test.int, scientific  = F)
  cat(paste0("\n Global integral difference test of (2000-2012) and (2013-2025) with multiplier bootstrap approximation:"))
  print(test.int)
  
  
  
  # not centered difference process
  test.delta  = readRDS(paste0("weather temperature/trend analysis/trend results/dense_quantile_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.delta  = p.adjust(p.val, method = "holm")
  test.delta  = format(test.delta, scientific  = F)
  cat(paste0("\n Global difference test of (2000-2012) and (2013-2025) with multiplier bootstrap approximation:"))
  print(test.delta)
  
  
  # centered difference process
  test.center.delta = readRDS(paste0("weather temperature/trend analysis/trend results/dense_P_quantile_",data.example[[k]],".rds"))
  p.val             = data.frame(t(test.center.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.center.delta = p.adjust(p.val, method = "holm")
  test.center.delta = format(test.center.delta, scientific  = F)
  cat(paste0("\n Global difference test of (1952-1962) and (1963-1972) with multiplier bootstrap approximation:"))
  print(test.center.delta)
  
  
  # Summary of all test: (2000-2012) with (2013-2025)
  years_2000_2012_with_2013_2025 = rbind(test.int,test.delta,test.center.delta)
  
  
  
  
  
  list.tests = list()
  
  list.tests$years_2000_2025_with_1952_1972 = years_2000_2025_with_1952_1972
  
  list.tests$robustness_15_check_years_2000_2025_with_1952_1972  = robustness_15_check_years_2000_2025_with_1952_1972
  list.tests$robustness_075_check_years_2000_2025_with_1952_1972 = robustness_075_check_years_2000_2025_with_1952_1972
  
  list.tests$years_1952_1962_with_1963_1972 = years_1952_1962_with_1963_1972
  list.tests$years_2000_2012_with_2013_2025 = years_2000_2012_with_2013_2025
  
  
  
  filename = paste0("weather temperature/global tests/global_tests_",data.example[[k]],".rds") 
  saveRDS(list.tests, filename)
  
  
}





#####################################################
### Change of station: Global test over 12 months ###
#####################################################


for(k in 1:2){
  
  
  data.example = list("Berlin", "Frankfurt_Main")      
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))
  station      = list("00433_00424","01425_01421") 
  
  #Integral: approximation with MB
  test.int = readRDS(paste0("weather temperature/trend analysis/change results/sparse_",station[[k]],"_quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.int    = p.adjust(p.val, method = "holm")
  test.int    = format(test.int, scientific  = F)
  if(k == 1){ periods = "(1952-1970) and (1971-1972)"}else{ periods = "(1952-1961) and (1962-1972)"}
  cat(paste0("\n Global integral difference test of ",periods," with multiplier bootstrap approximation:"))
  print(test.int)
  
  
  
  # not centered difference process
  test.delta  = readRDS(paste0("weather temperature/trend analysis/change results/sparse_",station[[k]],"_quantile_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.delta  = p.adjust(p.val, method = "holm")
  test.delta  = format(test.delta, scientific  = F)
  cat(paste0("\n Global difference test of ",periods," with multiplier bootstrap approximation:"))
  print(test.delta)
  
  
  # centered difference process
  test.centered.delta = readRDS(paste0("weather temperature/trend analysis/change results/sparse_",station[[k]],"_P_quantile_",data.example[[k]],".rds"))
  p.val = data.frame(t(test.centered.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.centered.delta = p.adjust(p.val, method = "holm")
  test.centered.delta= format(test.centered.delta, scientific  = F)
  cat(paste0("\n Global difference test of ",periods," with multiplier bootstrap approximation:"))
  print(test.centered.delta)

  
  list.tests = rbind(test.int,test.delta,test.centered.delta)
  
  
  filename = paste0("weather temperature/global tests/change_",station[[k]],"_global_tests_",data.example[[k]],".rds") 
  saveRDS(list.tests, filename)
  
}










#########################################################
### Relocation of station: Global test over 12 months ###
#########################################################


for(k in 1:3){
  
  
  data.example = list("Berlin", "Frankfurt_Main", "Hamburg") 
  station      = list("00433", "01420", "01975")
  regime       = "sparse"; if(k == 2) regime = "dense"
  periods      = list("(1952-1957) and (1957-1972)", "(2000-2014) and (2015-2025)", "(1952-1968) and (1969-1972)")
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))
  
  
  #Integral: approximation with MB
  test.int    = readRDS(paste0("weather temperature/trend analysis/relocation results/",regime,"_",station[[k]],"_quantile_integral_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.int$pval), row.names = c("H0: integral = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.int    = p.adjust(p.val, method = "holm")
  test.int    = format(test.int, scientific  = F)
  cat(paste0("\n Global integral difference test of ",periods[[k]]," with multiplier bootstrap approximation:"))
  print(test.int)
  
  
  
  # not centered difference process
  test.delta  = readRDS(paste0("weather temperature/trend analysis/relocation results/",regime,"_",station[[k]],"_quantile_",data.example[[k]],".rds"))
  p.val       = data.frame(t(test.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.delta  = p.adjust(p.val, method = "holm")
  test.delta  = format(test.delta, scientific  = F)
  cat(paste0("\n Global difference test of ",periods[[k]]," with multiplier bootstrap approximation:"))
  print(test.delta)
  
  
  
  # centered difference process
  test.centered.delta = readRDS(paste0("weather temperature/trend analysis/relocation results/",regime,"_",station[[k]],"_P_quantile_",data.example[[k]],".rds"))
  p.val               = data.frame(t(test.centered.delta[2,]), row.names = c("H0: delta = 0"))|> setNames(month.name)
  p.val[p.val == 0] = 1/10000
  test.centered.delta = p.adjust(p.val, method = "holm")
  test.centered.delta = format(test.centered.delta, scientific  = F)
  cat(paste0("\n Global difference test of ",periods[[k]]," with multiplier bootstrap approximation:"))
  print(test.centered.delta)
  
  
  list.tests = rbind(test.int,test.delta,test.centered.delta)
  
  
  filename = paste0("weather temperature/global tests/relocation_",station[[k]],"_global_tests_",data.example[[k]],".rds") 
  saveRDS(list.tests, filename)
  
}


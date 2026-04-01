library(MASS)
library(ggplot2)
library(reshape2)
library(locpol)
library(interp)
library(stats)
library(future.apply)
library(parallel)
library(biLocPol)
library(plotly)
library(tidyr)
library(goffda)
library(crayon)
library(future)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)


source("Mod_biLocPol.R")

source("Simulation/Simulation_Mod_biLocPol.R") #Adapted to grid from Berger and Holzmann (2025)
source("Simulation/function_data_generating.R")
source("Simulation/function_bandwidth_selection.R")



##########################################
# Change from H0 to H1           #########
#                                #########
# - TRUE:  delta constant        #########
# - FALSE: delta not constant    #########
##########################################
                          ################
delta_constant = TRUE     ################
                          ################
##########################################
if(delta_constant){       ################
  folder = "Delta = 0"    ################
}else{                    ################
  folder = "Delta = 0"    ################
}                         ################
##########################################


# Setting Parameters 
C_inv = 1.2
N  = c(15,25,50,75,100,150,200,250,300,350,400)
Nd = N * C_inv
P  = c(25,50,75)
pd = 100

#####################################################
# Computing optimal bandwidth for p = 25, 50 and 75 #
#####################################################


bw.list.p.25   = matrix(NA, nrow = 11, ncol = length(N),
                        dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.list.p.25) = N
bw.sample.list.p.25   = list()
bw.optim.list.p.25   = matrix(NA, nrow = 11, ncol = length(N),
                              dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.optim.list.p.25) = N

bw.list.p.50   = matrix(NA, nrow = 11, ncol = length(N),
                        dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.list.p.50) = N
bw.sample.list.p.50   = list()
bw.optim.list.p.50   = matrix(NA, nrow = 11, ncol = length(N),
                              dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.optim.list.p.50) = N

bw.list.p.75   = matrix(NA, nrow = 11, ncol = length(N),
                        dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.list.p.75) = N
bw.sample.list.p.75   = list()
bw.optim.list.p.75   = matrix(NA, nrow = 11, ncol = length(N),
                              dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.optim.list.p.75) = N

bw.list.pd.100   = matrix(NA, nrow = 11, ncol = length(Nd),
                          dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.list.pd.100) = Nd
bw.sample.list.pd.100   = list()
bw.optim.list.pd.100   = matrix(NA, nrow = 11, ncol = length(Nd),
                                dimnames = list(paste0("lag ", 0:10))) |> 
                                                              as.data.frame()
colnames(bw.optim.list.pd.100) = Nd


file25         = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_25_Bandwidth.rds")
file25.sample  = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_25_sample.rds")
file25.optim   = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_25_Bandwidth_optim.rds")

file50         = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_50_Bandwidth.rds")
file50.sample  = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_50_sample.rds")
file50.optim   = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_50_Bandwidth_optim.rds")

file75         = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_75_Bandwidth.rds")
file75.sample  = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_75_sample.rds")
file75.optim   = paste0("Simulation/optimal bandwidth/Results/lag_Gamma_s_75_Bandwidth_optim.rds")

file100        = paste0("Simulation/optimal bandwidth/Results/",folder,"/lag_Gamma_d_100_Bandwidth.rds")
file100.sample = paste0("Simulation/optimal bandwidth/Results/",folder,"/lag_Gamma_d_100_sample.rds")
file100.optim  = paste0("Simulation/optimal bandwidth/Results/",folder,"/lag_Gamma_d_100_Bandwidth_optim.rds")



# Optional: For parallelizing computations ###################
options(future.globals.maxSize = 64 * 1024^3)               ##  
plan(multisession, workers = 60) # MaRC3a: partition=normal ##
##############################################################

for(i in 2:length(N)){
  
  print(paste0("Start: n = ",N[i]))
  
  h.seq              = lapply(1, function(l){seq(1, 3/P[1], -0.05)})[[1]]
  bw.s.25.sample     = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[1], h.seq, K = 5, gap = 5, max.lag = 0, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.25            = mean(bw.s.25.sample)
  bw.s.lag.25.sample = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[1], h.seq, K = 5, gap = 5, max.lag = 4, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.lag.25        = apply(bw.s.lag.25.sample,1,mean)
  
  bw.list.p.25[,i]         = round(c(bw.s.25,bw.s.lag.25), digits = 2)
  bw.sample.list.p.25[[i]] = rbind(bw.s.25.sample,bw.s.lag.25.sample)
  
  h.seq              = lapply(1, function(l){seq(1, 4/P[1], -0.01)})[[1]]
  bw.s.25.optim      = h.optim.cov(N = 1000, n = N[i], p = P[1], p.eval = 75, h.seq, max.lag = 0,  m = 1, w.parallel = T, 
                                   func = mu)
  bw.s.lag.25.optim  = h.optim.cov(N = 1000, n = N[i], p = P[1], p.eval = 75, h.seq, max.lag = 4, m = 1, w.parallel = T, 
                                   func = mu) 
  
  bw.optim.list.p.25[,i]   = c(bw.s.25.optim,bw.s.lag.25.optim)
  
  saveRDS(round(bw.list.p.25, digits = 2),file25)
  saveRDS(bw.sample.list.p.25,file25.sample) 
  saveRDS(bw.optim.list.p.25,file25.optim)
  
  print(paste0("done: p = ",P[1]))
  
  
  
  h.seq              = lapply(1, function(l){seq(1, 3/P[2], -0.05)})[[1]]
  bw.s.50.sample     = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[2], h.seq, K = 5, gap = 5, max.lag = 0, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.50            = mean(bw.s.50.sample)
  bw.s.lag.50.sample = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[2], h.seq, K = 5, gap = 5, max.lag = 4, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.lag.50        = apply(bw.s.lag.50.sample,1,mean)
  
  bw.list.p.50[,i]         = c(bw.s.50,bw.s.lag.50)
  bw.sample.list.p.50[[i]] = rbind(bw.s.50.sample,bw.s.lag.50.sample)
  
  h.seq              = lapply(1, function(l){seq(1, 4.6/P[2], -0.01)})[[1]]
  bw.s.50.optim      = h.optim.cov(N = 1000, n = N[i], p = P[2], p.eval = 75, h.seq, max.lag = 0,  m = 1, w.parallel = T, 
                                   func = mu)
  bw.s.lag.50.optim  = h.optim.cov(N = 1000, n = N[i], p = P[2], p.eval = 75, h.seq, max.lag = 4, m = 1, w.parallel = T, 
                                   func = mu)
  
  bw.optim.list.p.50[,i]   = c(bw.s.50.optim,bw.s.lag.50.optim)
  
  saveRDS(round(bw.list.p.50, digits = 2),file50)
  saveRDS(bw.sample.list.p.50,file50.sample)
  saveRDS(bw.optim.list.p.50,file50.optim)
  
  print(paste0("done: p = ",P[2]))
  
  
  
  h.seq              = lapply(1, function(l){seq(1, 3/P[3], -0.05)})[[1]]
  bw.s.75.sample     = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[3], h.seq, K = 5, gap = 5, max.lag = 0, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.75            = mean(bw.s.75.sample)
  bw.s.lag.75.sample = k.fold.hc.cv.cov(N = 1000, n = N[i], p = P[3], h.seq, K = 5, gap = 5, max.lag = 4, m = 1, w.parallel = T, 
                                        func = mu)
  bw.s.lag.75        = apply(bw.s.lag.75.sample,1,mean)
  
  bw.list.p.75[,i]         = c(bw.s.75,bw.s.lag.75)
  bw.sample.list.p.75[[i]] = rbind(bw.s.75.sample,bw.s.lag.75.sample)
  
  h.seq              = lapply(1, function(l){seq(1, 4.6/P[3], -0.01)})[[1]]
  bw.s.75.optim      = h.optim.cov(N = 1000, n = N[i], p = P[3], p.eval = 75, h.seq, max.lag = 0,  m = 1, w.parallel = T, 
                                   func = mu)
  bw.s.lag.75.optim  = h.optim.cov(N = 1000, n = N[i], p = P[3], p.eval = 75, h.seq, max.lag = 4, m = 1, w.parallel = T, 
                                   func = mu)  
  
  bw.optim.list.p.75[,i]   = c(bw.s.75.optim,bw.s.lag.75.optim)
  
  saveRDS(round(bw.list.p.75, digits = 2),file75)
  saveRDS(bw.sample.list.p.75,file75.sample)
  saveRDS(bw.optim.list.p.75,file75.optim)
  
  print(paste0("done: p = ",P[3]))
  
  
  
  print(paste0("Start: nd = ",Nd[i]))
  
  h.seq              = lapply(1, function(l){seq(1, 3/pd, -0.05)})[[1]]
  bw.d.100.sample     = k.fold.hc.cv.cov(N = 1000, n = Nd[i], p = pd, h.seq, K = 5, gap = 5, max.lag = 0, m = 1, w.parallel = T, 
                                         func = function(x){mu_d(x,delta_constant)})
  bw.d.100            = mean(bw.d.100.sample)
  bw.d.lag.100.sample = k.fold.hc.cv.cov(N = 1000, n = Nd[i], p = pd, h.seq, K = 5, gap = 5, max.lag = 4, m = 1, w.parallel = T, 
                                         func = function(x){mu_d(x,delta_constant)})
  bw.d.lag.100        = apply(bw.d.lag.100.sample,1,mean)
  
  bw.list.pd.100[,i]         = c(bw.d.100,bw.d.lag.100)
  bw.sample.list.pd.100[[i]] = rbind(bw.d.100.sample,bw.d.lag.100.sample)
  
  h.seq              = lapply(1, function(l){seq(1, 3/pd, -0.05)})[[1]]
  bw.d.100.optim      = h.optim.cov(N = 1000, n = Nd[i], p = pd, p.eval = 75, h.seq, max.lag = 0, m = 1, w.parallel = T, 
                                    func = mu)
  bw.d.lag.100.optim  = h.optim.cov(N = 1000, n = Nd[i], p = pd, p.eval = 75, h.seq, max.lag = 4, m = 1, w.parallel = T, 
                                    func = mu)  
  
  bw.optim.list.pd.100[,i]   = c(bw.d.100.optim,bw.d.lag.100.optim)
  
  saveRDS(round(bw.list.pd.100, digits = 2),file100)
  saveRDS(bw.sample.list.pd.100,file100.sample)
  saveRDS(bw.optim.list.pd.100,file100.optim)
  
  print(paste0("done: pd = ",pd))
}

plan(sequential)










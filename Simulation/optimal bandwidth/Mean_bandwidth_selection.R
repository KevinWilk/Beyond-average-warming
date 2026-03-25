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

source("Simulation_Mod_biLocPol.R") #Adapted to grid from Berger and Holzmann (2025)
source("Simulation/function_data_generating.R")
source("Simulation/function_bandwidth_selection.R")



##########################################
# Change from H0 to H1           #########
#                                #########
# - TRUE:  delta constant        #########
# - FALSE: delta not constant    #########
##########################################
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

bw.optim.list = list()
bw.optim.list = matrix(NA, nrow = 8, ncol = length(N),
                       dimnames = list(c("dense.optim.100","dense.hv.100",
                                         "delta.optim.25","delta.hv.25",
                                         "delta.optim.50","delta.hv.50",
                                         "delta.optim.75","delta.hv.75"))) |> 
                                                                      as.data.frame()
colnames(bw.optim.list) = N



# Optional: For parallelizing computations ###############
options(future.globals.maxSize = 20 * 1024^3)           ##  
plan(multisession, workers = future::availableCores()-2)##
##########################################################

for(i in 1:length(N)){
  
  n  = N[i]
  nd = Nd[i]
  
  print(paste0("Start: n = ",N[i]))
  
  h.seq           = seq(0.08, 1, 0.01)
  bw.d.optim.list = h.optim(function(x){mu_d(x,delta_constant)},nd, pd, p.eval = 75, N = 1000, h.seq = h.seq, m = 2)
  bw.d.optim      = h.seq[which(bw.d.optim.list == min(bw.d.optim.list))]
  h.seq           = if(Nd[i] >= 60){seq(0.08, 1, 0.02)}else{seq(0.14, 1, 0.02)}
  bw.d.sample     = k.fold.hc.cv(N = 1000, n = nd, p = pd, K = 5 , gap = 5 , h.seq,  m = 2, 
                                 func = function(x){mu_d(x,delta_constant)})
  bw.d            = mean(bw.d.sample)
  
  print(paste0("done: dense for pd = 100"))
  
  h.seq                  = seq(0.08, 1, 0.01)
  bw.delta.optim.75.list = h.optim(function(x){delta(x,delta_constant)}, N[i], P[3], p.eval = 75, N = 1000, h.seq = h.seq, m = 2)
  bw.delta.optim.75      = h.seq[which(bw.delta.optim.75.list == min(bw.delta.optim.75.list))]
  h.seq                  = if(N[i] >= 50){seq(0.08, 1, 0.02)}else{seq(0.14, 1, 0.02)}
  bw.delta.sample.75     = k.fold.hc.cv(N = 1000, n = N[i], p = P[3], K = 5 , gap = 5 , h.seq,  m = 2, func = mu, diff = T,nd = nd, pd = pd, h.dense = bw.d, 
                                        func2 = function(x){mu_d(x,delta_constant)})
  bw.delta.75            = mean(bw.delta.sample.75)
  print(paste0("done: delta for p = ",P[3]))
  
  h.seq                  = seq(0.1, 1, 0.01)
  bw.delta.optim.50.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[2], p.eval = 75, N = 1000, h.seq = h.seq, m = 2)
  bw.delta.optim.50      = h.seq[which(bw.delta.optim.50.list == min(bw.delta.optim.50.list))]
  h.seq                  = if(N[i] >= 50){seq(0.1, 1, 0.02)}else{seq(0.16, 1, 0.02)}
  bw.delta.sample.50     = k.fold.hc.cv(N = 1000, n = N[i], p = P[2], K = 5 , gap = 5 , h.seq,  m = 2, func = mu, diff = T,nd = nd, pd = pd, h.dense = bw.d, 
                                        func2 = function(x){mu_d(x,delta_constant)})
  bw.delta.50            = mean(bw.delta.sample.50)
  print(paste0("done: delta for p = ",P[2]))
  
  h.seq                  = seq(0.1, 1, 0.01)
  bw.delta.optim.25.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[1], p.eval = 75, N = 1000, h.seq = h.seq, m = 2)
  bw.delta.optim.25      = h.seq[which(bw.delta.optim.25.list == min(bw.delta.optim.25.list))]
  h.seq                  = if(N[i] >= 50){seq(0.1, 1, 0.02)}else{seq(0.16, 1, 0.02)}
  bw.delta.sample.25     = k.fold.hc.cv(N = 1000, n = N[i], p = P[1], K = 5 , gap = 5 , h.seq,  m = 2, func = mu, diff = T,nd = nd, pd = pd, h.dense = bw.d, 
                                        func2 = function(x){mu_d(x,delta_constant)})
  bw.delta.25            = mean(bw.delta.sample.25)
  print(paste0("done: delta for p = ",P[1]))
  
  
  bw.optim.list[ , i] = round(c(bw.d.optim, bw.d, bw.delta.optim.25, bw.delta.25, bw.delta.optim.50, bw.delta.50, bw.delta.optim.75, bw.delta.75), digits = 2)
  
  
  file = paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds")
  saveRDS(bw.optim.list,file)
  
}


plan(sequential)




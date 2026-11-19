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
##########################################
                          ################
delta_constant = FALSE    ################
                          ################
##########################################
if(delta_constant){       ################
  folder = "Delta = 0"    ################
}else{                    ################
  folder = "Delta != 0"   ################
}                         ################
##########################################


# Setting Parameters 
C_inv = 1.2
N  = c(25,50,100,200,300,400,500,600,700,800,900,1000)
Nd = N * C_inv
P  = c(25,50,75,15,16,18,20)
pd = 100

#####################################################
# Computing optimal bandwidth for p = 25, 50 and 75 #
#####################################################

#bw.optim.list = list()
#bw.optim.list = matrix(NA, nrow = 30, ncol = length(N),
#                       dimnames = list(c("dense.optim.100","dense.hv.100",
                                         
#                                         "delta.optim.25" ,"delta.hv.25",
#                                         "delta.optim.50" ,"delta.hv.50",
#                                         "delta.optim.75" ,"delta.hv.75",
                                         
#                                         "sparse.optim.25" ,"sparse.hv.25",
#                                         "sparse.optim.50" ,"sparse.hv.50",
#                                         "sparse.optim.75" ,"sparse.hv.75",
                                         
#                                         "EXTRA: delta.optim.15" ,"EXTRA: delta.hv.15",
#                                         "EXTRA: sparse.optim.15" ,"EXTRA: sparse.hv.15",
                                         
#                                         "EXTRA: delta.optim.16" ,"EXTRA: delta.hv.16",
#                                         "EXTRA: sparse.optim.16" ,"EXTRA: sparse.hv.16",
                                         
#                                         "EXTRA: delta.optim.18" ,"EXTRA: delta.hv.18",
#                                         "EXTRA: sparse.optim.18" ,"EXTRA: sparse.hv.18",
                                         
#                                         "EXTRA: delta.optim.20" ,"EXTRA: delta.hv.20",
#                                         "EXTRA: sparse.optim.20" ,"EXTRA: sparse.hv.20"))) |> 
  
#                                                                      as.data.frame()
#colnames(bw.optim.list) = N


file = paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds")
bw.optim.list = readRDS(file)


# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 30 * 1024^3)                ##   
plan(multisession, workers = 30) # MaRC3a: partition=normal  ##
###############################################################

for(i in 1:length(N)){
  
  n  = N[i]
  nd = Nd[i]
  
  print(paste0("Start: n = ",N[i]))

  
  
  # Bandwidth selection of mu[d] by densely sampled sets
  h.seq           = seq(0.08, 1, 0.01)
  bw.d.optim.list = h.optim(function(x){mu_d(x,delta_constant)},nd, pd, p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.d.optim      = bw.d.optim.list$h.min
  h.seq           = if(Nd[i] >= 60){seq(0.08, 1, 0.01)}else{seq(0.14, 1, 0.01)}
  bw.d.sample     = cv.sim(N = 1000, n = nd, p = pd, h.seq, K = 5, gap = 5,  m = 2, R = 5, func = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.d            = mean(bw.d.sample)
  
  print(paste0("done: mu[d] for pd = 100"))
  
  
  # Bandwidth selection of delta by sparsely sampled as well as densely sampled sets
  h.seq                  = seq(0.08, 1, 0.01)
  bw.delta.optim.75.list = h.optim(function(x){delta(x,delta_constant)}, N[i], P[3], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.75      = bw.delta.optim.75.list$h.min
  h.seq                  = if(N[i] >= 50){seq(0.08, 1, 0.01)}else{seq(0.14, 1, 0.01)}
  bw.delta.sample.75     = cv.sim(N = 1000, n = N[i], p = P[3], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.75            = mean(bw.delta.sample.75)
  print(paste0("done: delta for p = ",P[3], ", h = ",bw.delta.75, ", sd = ", sd(bw.delta.sample.75) ))
  
  h.seq                  = seq(0.1, 1, 0.01)
  bw.delta.optim.50.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[2], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.50      = bw.delta.optim.50.list$h.min
  h.seq                  = if(N[i] >= 50){seq(0.08, 1, 0.01)}else{seq(0.16, 1, 0.01)}
  bw.delta.sample.50     = cv.sim(N = 1000, n = N[i], p = P[2], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.50            = mean(bw.delta.sample.50)
  print(paste0("done: delta for p = ",P[2], ", h = ",bw.delta.50, ", sd = ", sd(bw.delta.sample.50) ))

    
  h.seq                  = seq(0.1, 1, 0.01)
  bw.delta.optim.25.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[1], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.25      = bw.delta.optim.25.list$h.min
  h.seq                  = if(N[i] >= 50){seq(0.09, 1, 0.01)}else{seq(0.16, 1, 0.01)}
  bw.delta.sample.25     = cv.sim(N = 1000, n = N[i], p = P[1], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.25            = mean(bw.delta.sample.25)
  print(paste0("done: delta for p = ",P[1], ", h = ",bw.delta.25, ", sd = ", sd(bw.delta.sample.25) ))
  
  
    
    
    
  # Bandwidth selection of mu[s] only by sparse sample
  h.seq              = seq(0.08, 1, 0.01)
  bw.s.optim.75.list = h.optim(mu, N[i], P[3], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.75      = bw.s.optim.75.list$h.min
  h.seq              = if(N[i] >= 50){seq(0.08, 1, 0.01)}else{seq(0.14, 1, 0.01)}
  bw.s.sample.75     = cv.sim(N = 1000, n = N[i], p = P[3], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.75            = mean(bw.s.sample.75)
  print(paste0("done: mu[s] for p = ",P[3]))
  
  h.seq              = seq(0.1, 1, 0.01)
  bw.s.optim.50.list = h.optim(mu ,N[i], P[2], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.50      = bw.s.optim.50.list$h.min
  h.seq              = if(N[i] >= 50){seq(0.09, 1, 0.01)}else{seq(0.16, 1, 0.01)}
  bw.s.sample.50     = cv.sim(N = 1000, n = N[i], p = P[2], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.50            = mean(bw.s.sample.50)
  print(paste0("done: mu[s] for p = ",P[2]))
  
  h.seq              = seq(0.1, 1, 0.01)
  bw.s.optim.25.list = h.optim(mu ,N[i], P[1], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.25      = bw.s.optim.25.list$h.min
  h.seq              = if(N[i] >= 50){seq(0.09, 1, 0.01)}else{seq(0.16, 1, 0.01)}
  bw.s.sample.25     = cv.sim(N = 1000, n = N[i], p = P[1], h.seq, K = 5, gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.25            = mean(bw.s.sample.25)
  print(paste0("done: mu[s] for p = ",P[1]))


    
  # EXTRA #####################################################
  
  h.seq                  = if(N[i] >= 50){seq(0.18, 1, 0.01)}else{seq(0.21, 1, 0.01)}
  bw.delta.optim.15.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[4], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.15      = bw.delta.optim.15.list$h.min
  bw.delta.sample.15     = cv.sim(N = 1000, n = N[i], p = P[4], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, diff = T,nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.15            = mean(bw.delta.sample.15)
  print(paste0("done: delta for p = ",P[4]))
  
  h.seq              = if(N[i] >= 50){seq(0.18, 1, 0.01)}else{seq(0.2, 1, 0.01)}
  bw.s.optim.15.list = h.optim(mu ,N[i], P[4], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.15      = bw.s.optim.15.list$h.min
  bw.s.sample.15     = cv.sim(N = 1000, n = N[i], p = P[4], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.15            = mean(bw.s.sample.15)
  print(paste0("done: mu[s] for p = ",P[4]))
  
  
  h.seq                  = if(N[i] >= 50){seq(0.17, 1, 0.01)}else{seq(0.2, 1, 0.01)}
  bw.delta.optim.16.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[5], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.16      = bw.delta.optim.16.list$h.min
  bw.delta.sample.16     = cv.sim(N = 1000, n = N[i], p = P[5], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.16            = mean(bw.delta.sample.16)
  print(paste0("done: delta for p = ",P[5]))
  
  h.seq              = if(N[i] >= 50){seq(0.17, 1, 0.01)}else{seq(0.2, 1, 0.01)}
  bw.s.optim.16.list = h.optim(mu ,N[i], P[5], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.16      = bw.s.optim.16.list$h.min
  bw.s.sample.16     = cv.sim(N = 1000, n = N[i], p = P[5], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.16            = mean(bw.s.sample.16)
  print(paste0("done: mu[s] for p = ",P[5]))
  

  h.seq              = if(N[i] >= 50){seq(0.16, 1, 0.01)}else{seq(0.18, 1, 0.01)}
  bw.delta.optim.18.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[6], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.18      = bw.delta.optim.18.list$h.min
  bw.delta.sample.18     = cv.sim(N = 1000, n = N[i], p = P[6], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.18            = median(bw.delta.sample.18)
  print(paste0("done: delta for p = ",P[6]))
  
  h.seq              = if(N[i] >= 50){seq(0.16, 1, 0.01)}else{seq(0.18, 1, 0.01)}
  bw.s.optim.18.list = h.optim(mu ,N[i], P[6], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.18      = bw.s.optim.18.list$h.min
  bw.s.sample.18     = cv.sim(N = 1000, n = N[i], p = P[6], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.18            = mean(bw.s.sample.18)
  print(paste0("done: mu[s] for p = ",P[6]))
  

  h.seq              = if(N[i] >= 50){seq(0.15, 1, 0.01)}else{seq(0.12, 1, 0.01)}
  bw.delta.optim.20.list = h.optim(function(x){delta(x,delta_constant)} ,N[i], P[7], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.delta.optim.20      = bw.delta.optim.20.list$h.min
  bw.delta.sample.20     = cv.sim(N = 1000, n = N[i], p = P[7], h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, diff = T, nd = nd, pd = pd, h.dense = bw.d, func2 = function(x){mu_d(x,delta_constant)}, norm = "L2-norm")
  bw.delta.20            = mean(bw.delta.sample.20)
  print(paste0("done: delta for p = ",P[7]))
  
  h.seq              = if(N[i] >= 50){seq(0.15, 1, 0.01)}else{seq(0.12, 1, 0.01)}
  bw.s.optim.20.list = h.optim(mu ,N[i], P[7], p.eval = 75, N = 1000, h.seq = h.seq, m = 2, norm = "L2-norm")
  bw.s.optim.20      = bw.s.optim.20.list$h.min
  bw.s.sample.20     = cv.sim(N = 1000, n = N[i], p = P[7],  h.seq, K = 5 , gap = 5,  m = 2, R = 5, func = mu, norm = "L2-norm")
  bw.s.20            = mean(bw.s.sample.20)
  print(paste0("done: mu[s] for p = ",P[7]))
    
  #################################################################
  
  bw.optim.list[ , i] = round(c(bw.d.optim, bw.d, 
                                
                                bw.delta.optim.25, bw.delta.25, 
                                bw.delta.optim.50, bw.delta.50, 
                                bw.delta.optim.75, bw.delta.75,
                                
                                bw.s.optim.25, bw.s.25, 
                                bw.s.optim.50, bw.s.50, 
                                bw.s.optim.75, bw.s.75,
                                
                                bw.delta.optim.15, bw.delta.15, 
                                bw.s.optim.15, bw.s.15,
                                
                                bw.delta.optim.16, bw.delta.16, 
                                bw.s.optim.16, bw.s.16,
                                
                                bw.delta.optim.18, bw.delta.18, 
                                bw.s.optim.18, bw.s.18,
                                
                                bw.delta.optim.20, bw.delta.20, 
                                bw.s.optim.20, bw.s.20), 
                              
                              digits = 2)
  
  
  file = paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds")
  saveRDS(bw.optim.list,file)
  
}





plan(sequential)



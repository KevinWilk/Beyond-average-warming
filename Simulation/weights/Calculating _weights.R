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



# Setting Parameters 
C_inv = 1.2
N  = c(15,25,50,75,100,150,200,250,300,350,400)
Nd = N * C_inv
P  = c(25,50,75)
pd = 100
p.eval = 75



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



#######################################################
# Calculating weigths of bivariate covariance kernel: #
# for all b = -3,...,3                                #
#######################################################



# Loading the selected bandwidths in sparse setting
lag_Gamma_s_25_Bandwidth  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_25_Bandwidth.rds")
lag_Gamma_s_50_Bandwidth  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_50_Bandwidth.rds")
lag_Gamma_s_75_Bandwidth  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_75_Bandwidth.rds")


#################
# sparse sample #
#################


# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 64 * 1024^3)                ##   
plan(multisession, workers = 120) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################
  
for(j in 2:length(N)){

  for(i in 1:3){
    file.s = paste0("Simulation/weights/p ",P[i],"/full_w_lag0_",N[j],".rds")
    bandwidth = get(paste0("lag_Gamma_s_",P[i],"_Bandwidth"))
    w  = local.polynomial.weights(P[i], bandwidth[1,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "less", eval.type = "full", parallel.environment = F)
    saveRDS(w, file = file.s)
  
    file.s = paste0("Simulation/weights/p ",P[i],"/full_w_lag1_",N[j],".rds")
    bandwidth = get(paste0("lag_Gamma_s_",P[i],"_Bandwidth"))
    w  = local.polynomial.weights(P[i], bandwidth[2,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)
    saveRDS(w, file = file.s)
  
    file.s = paste0("Simulation/weights/p ",P[i],"/full_w_lag2_",N[j],".rds")
    bandwidth = get(paste0("lag_Gamma_s_",P[i],"_Bandwidth"))
    w  = local.polynomial.weights(P[i], bandwidth[3,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)
    saveRDS(w, file = file.s)
  
    file.s = paste0("Simulation/weights/p ",P[i],"/full_w_lag3_",N[j],".rds")
    bandwidth = get(paste0("lag_Gamma_s_",P[i],"_Bandwidth"))
    w  = local.polynomial.weights(P[i], bandwidth[4,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)
    saveRDS(w, file = file.s)
    rm(w,bandwidth,file.s); gc()
  }
}

plan(sequential)




################
# dense sample #
################
# - for pd 100 constant TRUE  set 'delta constant = TRUE'
# - for pd 100 constant FALSE set 'delta constant = FALSE'


# Loading the selected bandwidths in dense setting
lag_Gamma_d_100_Bandwidth = readRDS(paste0("Simulation/optimal bandwidth/Results/",if(delta_constant){"Delta = 0"}else{"Delta != 0"},"/lag_Gamma_d_100_Bandwidth.rds"))

# Optional: For parallelizing computations ###############
options(future.globals.maxSize = 20 * 1024^3)           ##  
plan(multisession, workers = future::availableCores()-2)##
##########################################################

for(j in 2:length(N)){

  file.d = paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag0_",N[j],".rds")
  bandwidth = get(paste0("lag_Gamma_d_",pd,"_Bandwidth"))
  wd = local.polynomial.weights(pd, bandwidth[1,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "less", eval.type = "full", parallel.environment = F)    
  saveRDS(wd, file = file.d)

  file.d = paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag1_",N[j],".rds")
  bandwidth = get(paste0("lag_Gamma_d_",pd,"_Bandwidth"))
  wd = local.polynomial.weights(pd, bandwidth[2,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)    
  saveRDS(wd, file = file.d)

  file.d = paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag2_",N[j],".rds")
  bandwidth = get(paste0("lag_Gamma_d_",pd,"_Bandwidth"))
  wd = local.polynomial.weights(pd, bandwidth[3,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)    
  saveRDS(wd, file = file.d)

  file.d = paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag3_",N[j],".rds")
  bandwidth = get(paste0("lag_Gamma_d_",pd,"_Bandwidth"))
  wd = local.polynomial.weights(pd, bandwidth[4,j], p.eval = p.eval, parallel = T, m = 1, del = 0, grid.type = "lesseq", eval.type = "full", parallel.environment = F)    
  saveRDS(wd, file = file.d) 
    
  rm(wd,bandwidth,file.d); gc()

}
    
plan(sequential)





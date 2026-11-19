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
source("Simulation/function_evaluation.R")



# Setting Parameters 
C_inv = 1.2
N  = c(25,50,100,200,300,400,500,600,700,800,900,1000)
Nd = N * C_inv
P  = c(25,50,75,15,16,18,20)
pd = 100
p.eval = 75

Repitition = 1000  # Number of data simulation repetitions 


##########################################
# Change from H0 to H1           #########
#                                #########
# - TRUE:  delta constant        #########
# - FALSE: delta not constant    #########
##########################################
##########################################
#                          ###############
#delta_constant = TRUE     ###############
#integral       = TRUE     ###############
#                          ###############
##########################################
##########################################
if(delta_constant){                  #####
  if(integral){                      #####
    folder = "H0 (centered)"         #####
  }else{folder = "H0 (not centered)"}#####
    ######################################
}else{                               #####
  if(integral){                      #####
    folder = "H1 (centered)"         #####
  }else{folder = "H1 (not centered)"}#####
}                                    #####
##########################################

if(scenario3){
    if(integral){
        folder = "Scenario 3 (centered)"
    }
    mu = function(x){
      x = 2*x - 1
      return(3*sin(1.5*pi*x)*exp(-2*abs(x)))
    }

    delta = function(x, constant = delta_constant){
      return(2*mu(x))
    }
}

sup.error.res = c()
sup.error.sep = c()

# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 30 * 1024^3)                ##   
plan(multisession, workers = 30) # MaRC3a: partition=mqtest ##
#plan(multisession, workers = 60) # MaRC3a: partition=normal ##
###############################################################

for(rep in 1:Repitition){
  
  
  
  #Gernate data sets
  x.s.design = (1:P[i] - 0.5)/P[i]
  Y.s       = t(mu(x.s.design) + 
                  t(sim.d.OU(N[j], t = x.s.design, rho_B = 0.5,tau = 0,sigma = 4)) + 
                    matrix(rnorm(length(x.s.design)*N[j], 0, 0.1),length(x.s.design), N[j]))
  
  xd.design = (1:pd - 0.5)/pd
  Y.d       = t(mu_d(xd.design,delta_constant) + 
                  t(sim.d.OU(Nd[j], t = xd.design, rho_B = 0.5,tau = 0,sigma = 4)) + 
                    matrix(rnorm(length(xd.design)*Nd[j], 0, 0.1),length(xd.design), Nd[j]))
  

  
  ##################################
  # Mean and Difference Estimation #
  ##################################
  
  # bandwidths from 5-Fold hv-block cross validation 
  h.d.seq  = seq(0.04, 1, 0.01)
  bw.d     = kfold.cv(Y.d, h.d.seq, K = 5, gap = 5, m = 2, R = 5, parallel = T, norm = "L2-norm")
  w.d      = locPolWeights(x = (1:ncol(Y.d) - 0.5)/ncol(Y.d), bw = bw.d, deg = 2, xeval= (1:ncol(Y.s) - 0.5)/ncol(Y.s), kernel=EpaK)$locWeig
  residual = t(t(Y.s) - c(w.d %*% apply(Y.d,2,mean)))
  h.s.seq  = seq(max(0.06,3/P[i]), 1, 0.01)
  
  #residuals based estimator
  bw.s    = kfold.cv(residual, h.s.seq, K = 5, gap = 5, m = 2, R = 5, parallel = T, norm = "L2-norm")
  est.cv  = est.results(Y.s,Y.d, c(bw.d, bw.s), p.eval = p.eval)
  est.res = if(integral){est.cv$delta$ESTIMATE - est.cv$delta_int$ESTIMATE }else{ est.cv$delta$ESTIMATE }

  # estimator by separately estimating both mean functions
  bw.sep      = kfold.cv(Y.s, h.s.seq, K = 5, gap = 5, m = 2, R = 5, parallel = T, norm = "L2-norm")
  est.cv.sep  = est.results(Y.s,Y.d, c(bw.d,bw.sep), p.eval = p.eval)  
  est.sep     = if(integral){ est.cv.sep$delta$ESTIMATE - est.cv.sep$delta_int$ESTIMATE  }else{ est.cv.sep$delta$ESTIMATE }
    
    
  ###################
  # sup-norm error  #
  ###################
    
  # true function 
  true.val      = delta((1:p.eval-0.5)/p.eval, constant = delta_constant)                                                                       
  res.error     = max(abs(true.val - est.cv$delta$ESTIMATE))
  sup.error.res = c(sup.error.res, res.error)
  sep.error     = max(abs(true.val - est.cv.sep$delta$ESTIMATE))
  sup.error.sep = c(sup.error.sep, sep.error)
  
  
  if(rep %% 10 == 0){print(rep)}
    

}

plan(sequential)


###################
### Saving data ### 
###################

file  = paste0("Simulation/Comparison/", folder, "/p ", P[i], "/dep_", N[j], "_rep_", Repitition, "_result.rds")

results = list( sup.error  = list(seperat.est  = sup.error.sep,
                                  residual.est =sup.error.res)
)

saveRDS(results,  file)















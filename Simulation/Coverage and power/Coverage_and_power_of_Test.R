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
N  = c(15,25,50,75,100,150,200,250,300,350,400)
Nd = N * C_inv
P  = c(25,50,75)
pd = 100
p.eval = 75

B          = 1000  # Number of repetitions of Multiplier Bootstrap
Repitition = 1000  # Number of data simulation repetitions 



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
##########################################
if(delta_constant){       ################
  folder = "Under H0"     ################
}else{                    ################
  folder = "Under H1"     ################
}                         ################
##########################################



# saving 1000 times calculated test statistic
emp_sample   = c()

# saving 1000 times estimated quantiles (from dependent Multiplier Bootstrap)
dep_q90_list     = c()
dep_q95_list     = c()
dep_q99_list     = c()

# Calculate incorrect rejection under the null hypothesis
dep.m90 = 0
dep.m95 = 0
dep.m99 = 0

if(!delta_constant){
  
  # correct rejection under the alternative
  
  dep.power90 = 0
  dep.power95 = 0
  dep.power99 = 0
}



# saving 1000 times estimated quantiles (from independent Multiplier Bootstrap)
ind_q90_list     = c()
ind_q95_list     = c()
ind_q99_list     = c()

# Calculate incorrect rejection under the null hypothesis
ind.m90 = 0
ind.m95 = 0
ind.m99 = 0

if(!delta_constant){
  
  # correct rejection under the alternative
  
  ind.power90 = 0
  ind.power95 = 0
  ind.power99 = 0
}


# Optional: For parallelizing computations ####################
options(future.globals.maxSize = 64 * 1024^3)                ##   
plan(multisession, workers = 120) # MaRC3a: partition=mqtest ##
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
  
  
  
  
  
  #########################################
  # Cumputing long run covariance kernel  #
  #########################################
  
  #sparse 
  w       = readRDS(paste0("Simulation/weights/p ",P[i],"/full_w_lag0_",N[j],".rds"))
  Gamma.s = eval.weights(w,observation.transformation(Y.s,grid.type = "less"))
  
  
  Gamma.s.part2 = Reduce(`+`,future_lapply(1:3, function(k) {
    
    w.lag       = readRDS(paste0("Simulation/weights/p ",P[i],"/full_w_lag",k,"_",N[j],".rds"))
    lag.Gamma.s = eval.weights(w.lag, observation.transformation(Y.s, lag = k, grid.type = "lesseq"), lag = k)
    lag.Gamma.s = (1 - k/(3 + 1) ) * (lag.Gamma.s + t(lag.Gamma.s))
    
    lag.Gamma.s 
  }))
  
  
  Gamma.s   = Gamma.s + Gamma.s.part2
  P.Gamma.s = diag(Gamma.s) - rowMeans(Gamma.s) - colMeans(Gamma.s) + mean(Gamma.s) 
  
  
  #dense 
  wd      = readRDS(paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag0_",N[j],".rds"))
  Gamma.d = eval.weights(wd,observation.transformation(Y.d,grid.type = "less"))
  
  
  Gamma.d.part2 = Reduce(`+`,future_lapply(1:3, function(k) {
    
    wd.lag      = readRDS(paste0("Simulation/weights/pd 100 constant ",delta_constant,"/full_w_lag",k,"_",N[j],".rds"))
    lag.Gamma.d = eval.weights(wd.lag, observation.transformation(Y.d, lag = k, grid.type = "lesseq"), lag = k)
    lag.Gamma.d = (1 - k/(3 + 1) ) * (lag.Gamma.d + t(lag.Gamma.d))
    
    lag.Gamma.d 
  }))
  
  
  Gamma.d   = Gamma.d + Gamma.d.part2
  P.Gamma.d = diag(Gamma.d) - rowMeans(Gamma.d) - colMeans(Gamma.d) + mean(Gamma.d) 
  
  
  ### long run kernel ################################# 
  Gamma   = diag(Gamma.s) + N[j]/Nd[j]* diag(Gamma.d) #
  P.Gamma = P.Gamma.s     + N[j]/Nd[j]* P.Gamma.d     # 
  #####################################################
  
  
  
  # Mean and Difference Estimation #
  bandwidth  = readRDS(paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds"))
  estimation = est.results(Y.s,Y.d,bandwidth[c(2,i*2+2),j], p.eval = p.eval)  
  
  
  delta.centered     = delta((1:p.eval-0.5)/p.eval, constant = delta_constant)-integrate(function(x){delta(x,constant = delta_constant)},0,1)$value 
  delta.centered.est = estimation$delta$ESTIMATE - estimation$delta_int$ESTIMATE
  
  
  ##################################
  # Dependent Multiplier Bootstrap #
  ##################################
  
  dep.sim = q.MB(Y.s, Y.d, estimation$sparse, estimation$dense, P.Gamma, (1:p.eval)/p.eval, bandwidth[c(2,i*2+2),j], B = B, depend = T, int = T, H0 = 0) 

  
  dep.q90 = quantile(dep.sim$sample, probs = 0.9,  Type = 2, na.rm = T)
  dep.q95 = quantile(dep.sim$sample, probs = 0.95, Type = 2, na.rm = T)
  dep.q99 = quantile(dep.sim$sample, probs = 0.99, Type = 2, na.rm = T)
  
  q90_list = c(dep_q90_list,dep.q90)
  q95_list = c(dep_q95_list,dep.q95)
  q99_list = c(dep_q99_list,dep.q99)

  
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*dep.q90) ){dep.m90 = dep.m90+1}
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*dep.q95) ){dep.m95 = dep.m95+1}
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*dep.q99) ){dep.m99 = dep.m99+1}
  
  if(!delta_constant){
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*dep.q90) ){dep.power90 = dep.power90+1}
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*dep.q95) ){dep.power95 = dep.power95+1}
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*dep.q99) ){dep.power99 = dep.power99+1}
  }
  
  
  
  
  ####################################
  # Independent Multiplier Bootstrap #
  ####################################
  
  ind.sim = q.MB(Y.s, Y.d, estimation$sparse, estimation$dense, P.Gamma, (1:p.eval)/p.eval, bandwidth[c(2,i*2+2),j], B = B, depend = F, int = T, H0 = 0)
  
  
  ind.q90 = quantile(ind.sim$sample, probs = 0.9,  Type = 2, na.rm = T)
  ind.q95 = quantile(ind.sim$sample, probs = 0.95, Type = 2, na.rm = T)
  ind.q99 = quantile(ind.sim$sample, probs = 0.99, Type = 2, na.rm = T)
  
  q90_list = c(ind_q90_list,ind.q90)
  q95_list = c(ind_q95_list,ind.q95)
  q99_list = c(ind_q99_list,ind.q99)
  
  
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*ind.q90) ){ind.m90 = ind.m90+1}
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*ind.q95) ){ind.m95 = ind.m95+1}
  if( any(abs(delta.centered.est - delta.centered) > sqrt(P.Gamma/(N[j]))*ind.q99) ){ind.m99 = ind.m99+1}
  
  if(!delta_constant){
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*ind.q90) ){ind.power90 = ind.power90+1}
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*ind.q95) ){ind.power95 = ind.power95+1}
    if( any(abs(delta.centered.est - 0) > sqrt(P.Gamma/(N[j]))*ind.q99) ){ind.power99 = ind.power99+1}
  }
  
  
  
  
  if(rep %% 10 == 0){print(rep)}
  
  emp_sample = c(emp_sample,sqrt(N[j])*max(abs(delta.centered.est/sqrt(P.Gamma))))
  
}

plan(sequential)





###################
### Saving data ### 
###################

# Dependent Multiplier Bootstrap

file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_coverage90.rds")
coverage      = (Repitition-dep.m90)/Repitition
saveRDS(coverage,file_coverage)
file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_coverage95.rds")
coverage      = (Repitition-dep.m95)/Repitition
saveRDS(coverage,file_coverage)
file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_coverage99.rds")
coverage      = (Repitition-dep.m99)/Repitition
saveRDS(coverage,file_coverage)

if(!delta_constant){
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_power90.rds")
  power      = dep.power90/Repitition
  saveRDS(power,file_power)
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_power95.rds")
  power      = dep.power95/Repitition
  saveRDS(power,file_power)
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_power99.rds")
  power      = dep.power99/Repitition
  saveRDS(power,file_power)
}

file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_q90_list.rds")
saveRDS(dep_q90_list,file_dep)
file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_q95_list.rds")
saveRDS(dep_q95_list,file_dep)
file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Repitition,"_q99_list.rds")
saveRDS(dep_q99_list,file_dep)




# Independent Multiplier Bootstrap

file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_coverage90.rds")
coverage      = (Repitition-ind.m90)/Repitition
saveRDS(coverage,file_coverage)
file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_coverage95.rds")
coverage      = (Repitition-ind.m95)/Repitition
saveRDS(coverage,file_coverage)
file_coverage = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_coverage99.rds")
coverage      = (Repitition-ind.m99)/Repitition
saveRDS(coverage,file_coverage)

if(!delta_constant){
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_power90.rds")
  power      = ind.power90/Repitition
  saveRDS(power,file_power)
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_power95.rds")
  power      = ind.power95/Repitition
  saveRDS(power,file_power)
  file_power = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_power99.rds")
  power      = ind.power99/Repitition
  saveRDS(power,file_power)
}

file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_q90_list.rds")
saveRDS(ind_q90_list,file_dep)
file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_q95_list.rds")
saveRDS(ind_q95_list,file_dep)
file_dep = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Repitition,"_q99_list.rds")
saveRDS(ind_q99_list,file_dep)





# Saving empirical quantile of test statistic from 1000 samples

q_emp      = quantile(sort(emp_sample),probs = 0.9, Type = 2, na.rm = T)
file_emp   = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/",N[j],"_rep_",Repitition,"_emp90.rds")
saveRDS(q_emp,file_emp)
q_emp      = quantile(sort(emp_sample),probs = 0.95, Type = 2, na.rm = T)
file_emp   = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/",N[j],"_rep_",Repitition,"_emp95.rds")
saveRDS(q_emp,file_emp)
q_emp      = quantile(sort(emp_sample),probs = 0.99, Type = 2, na.rm = T)
file_emp   = paste0("Simulation/coverage and power/",folder,"/p ",P[i],"/",N[j],"_rep_",Repitition,"_emp99.rds")
saveRDS(q_emp,file_emp)




rm(list = ls())








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
#delta_constant = TRUE     ###############
#integral       = TRUE     ###############
                          ################
##########################################
##########################################
##########################################
if(delta_constant){                  #####
  if(integral){                      #####
        folder = "H0 (centered)"     #####
  }else{folder = "H0 (not centered)"}#####
    ######################################
}else{                               #####
  if(integral){                      #####
        folder = "H1 (centered)"     #####
  }else{folder = "H1 (not centered)"}#####
}                                    #####
##########################################



# saving 1000 times calculated test statistics
emp.stat.L2 = c(); emp.stat.fix = c(); emp.stat.scaleL2 = c(); emp.stat.fix08 = c()

# saving 1000 times calculated CV bandwidths
h.s.list.L2 = c()
h.d.list.L2 = c()

# saving 1000 times estimated quantiles (from dependent Multiplier Bootstrap)
q90.L2.list = c(); q90.fix.list = c(); q90.scaleL2.list = c()
q95.L2.list = c(); q95.fix.list = c(); q95.scaleL2.list = c()
q99.L2.list = c(); q99.fix.list = c(); q99.scaleL2.list = c()

# counting separately each sample repetition
rej90.L2.list = c(); rej90.fix.list = c(); rej90.scaleL2.list = c()
rej95.L2.list = c(); rej95.fix.list = c(); rej95.scaleL2.list = c()
rej99.L2.list = c(); rej99.fix.list = c(); rej99.scaleL2.list = c()


# Calculate incorrect rejection under the null hypothesis
m90.L2 = 0; m90.fix = 0; m90.scaleL2 = 0
m95.L2 = 0; m95.fix = 0; m95.scaleL2 = 0
m99.L2 = 0; m99.fix = 0; m99.scaleL2 = 0

if(!delta_constant){
  
  # correct rejection under the alternative
  
  power90.L2 = 0; power90.fix = 0; power90.scaleL2 = 0
  power95.L2 = 0; power95.fix = 0; power95.scaleL2 = 0
  power99.L2 = 0; power99.fix = 0; power99.scaleL2 = 0
  
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
  
  kernel.diag = if(integral){ P.Gamma }else{ Gamma }
  

  
  
  # true function ##################################################################################################################
  true.val = if(integral){                                                                                                         #
    delta((1:p.eval-0.5)/p.eval, constant = delta_constant) - integrate(function(x){delta(x,constant = delta_constant)},0,1)$value #
  }else{                                                                                                                           #
    delta((1:p.eval-0.5)/p.eval, constant = delta_constant)}                                                                       #
  ##################################################################################################################################
  
  
  
  ##################################  
  # Mean and Difference Estimation #
  ################################## 
  
  
  # bandwidths from 5-Fold hv-block cross validation with L2-norm
  h.d.seq  = seq(0.04, 1, 0.01)
  bw.d.L2  = kfold.cv(Y.d, h.d.seq, K = 5, gap = 5, m = 2, R = 5, parallel = T, norm = "L2-norm")
  w.d      = locPolWeights(x = (1:ncol(Y.d) - 0.5)/ncol(Y.d), bw = bw.d.L2, deg = 2, xeval= (1:ncol(Y.s) - 0.5)/ncol(Y.s), kernel=EpaK)$locWeig
  residual = t(t(Y.s) - c(w.d %*% apply(Y.d,2,mean)))
  h.s.seq  = seq(max(0.06,3/P[i]), 1, 0.01)
  bw.s.L2  = kfold.cv(residual, h.s.seq, K = 5, gap = 5, m = 2, R = 5, parallel = T, norm = "L2-norm")
  est.L2  = est.results(Y.s,Y.d, c(bw.d.L2, bw.s.L2), p.eval = p.eval)
  obs.L2  = if(integral){est.L2$delta$ESTIMATE - est.L2$delta_int$ESTIMATE }else{ est.L2$delta$ESTIMATE }
  
  h.s.list.L2 = c(h.s.list.L2, bw.s.L2)
  h.d.list.L2 = c(h.d.list.L2, bw.d.L2)
  
  
  # Reference fix L2-norm bandwidth
  bandwidth = readRDS(paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds"))
  bw.d.fix  = bandwidth[2,j]
  bw.s.fix  = bandwidth[2+i*2,j]
  est.fix   = est.results(Y.s,Y.d, c(bw.d.fix, bw.s.fix), p.eval = p.eval)
  obs.fix   = if(integral){est.fix$delta$ESTIMATE - est.fix$delta_int$ESTIMATE }else{ est.fix$delta$ESTIMATE }

  
  # Reference CV L2-norm bandwidth * log(n^(-0.2))
  bw.d.scaleL2  = max(0.04,             round(bw.d.L2 * log(N[j])^(-0.2), digits = 2))
  bw.s.scaleL2  = max(max(0.06,3/P[i]), round(bw.s.L2 * log(N[j])^(-0.2), digits = 2))
  est.scaleL2   = est.results(Y.s,Y.d, c(bw.d.scaleL2, bw.s.scaleL2), p.eval = p.eval)
  obs.scaleL2   = if(integral){est.scaleL2$delta$ESTIMATE - est.scaleL2$delta_int$ESTIMATE }else{ est.scaleL2$delta$ESTIMATE }
  
  
  ##################################
  # Dependent Multiplier Bootstrap #
  ##################################
  
  sim.MB.L2      = q.MB(Y.s, Y.d, est.L2$sparse,      est.L2$dense,      kernel.diag, (1:p.eval - 0.5)/p.eval, c(bw.d.L2,bw.s.L2),           B = B, depend = T, int = integral) 
  sim.MB.fix     = q.MB(Y.s, Y.d, est.fix$sparse,     est.fix$dense,     kernel.diag, (1:p.eval - 0.5)/p.eval, c(bw.d.fix,bw.s.fix),         B = B, depend = T, int = integral)
  sim.MB.scaleL2 = q.MB(Y.s, Y.d, est.scaleL2$sparse, est.scaleL2$dense, kernel.diag, (1:p.eval - 0.5)/p.eval, c(bw.d.scaleL2,bw.s.scaleL2), B = B, depend = T, int = integral)
  
  
  q90.L2 = quantile(sim.MB.L2$sample, probs = 0.9,   type = 2, na.rm = T)
  q95.L2 = quantile(sim.MB.L2$sample, probs = 0.95,  type = 2, na.rm = T)
  q99.L2 = quantile(sim.MB.L2$sample, probs = 0.99,  type = 2, na.rm = T)
  
  q90.L2.list = c(q90.L2.list,q90.L2)
  q95.L2.list = c(q95.L2.list,q95.L2)
  q99.L2.list = c(q99.L2.list,q99.L2)
  
  
  
  q90.fix = quantile(sim.MB.fix$sample, probs = 0.9,   type = 2, na.rm = T)
  q95.fix = quantile(sim.MB.fix$sample, probs = 0.95,  type = 2, na.rm = T)
  q99.fix = quantile(sim.MB.fix$sample, probs = 0.99,  type = 2, na.rm = T)
  
  q90.fix.list = c(q90.fix.list,q90.fix)
  q95.fix.list = c(q95.fix.list,q95.fix)
  q99.fix.list = c(q99.fix.list,q99.fix)
  
  
  
  
  q90.scaleL2 = quantile(sim.MB.scaleL2$sample, probs = 0.9,   type = 2, na.rm = T)
  q95.scaleL2 = quantile(sim.MB.scaleL2$sample, probs = 0.95,  type = 2, na.rm = T)
  q99.scaleL2 = quantile(sim.MB.scaleL2$sample, probs = 0.99,  type = 2, na.rm = T)
  
  q90.scaleL2.list = c(q90.scaleL2.list,q90.scaleL2)
  q95.scaleL2.list = c(q95.scaleL2.list,q95.scaleL2)
  q99.scaleL2.list = c(q99.scaleL2.list,q99.scaleL2)
  
  
  
  # Width of confindence bands
  w.90.L2 = sqrt(kernel.diag / N[j]) * q90.L2
  w.95.L2 = sqrt(kernel.diag / N[j]) * q95.L2
  w.99.L2 = sqrt(kernel.diag / N[j]) * q99.L2
  
  w.90.fix = sqrt(kernel.diag / N[j]) * q90.fix
  w.95.fix = sqrt(kernel.diag / N[j]) * q95.fix
  w.99.fix = sqrt(kernel.diag / N[j]) * q99.fix
  
  w.90.scaleL2 = sqrt(kernel.diag / N[j]) * q90.scaleL2
  w.95.scaleL2 = sqrt(kernel.diag / N[j]) * q95.scaleL2
  w.99.scaleL2 = sqrt(kernel.diag / N[j]) * q99.scaleL2
  
  
  # L2-norm bandwidth
  rej90.L2 = 0; rej95.L2 = 0; rej99.L2 = 0
  if( any(abs(obs.L2 - true.val) > w.90.L2) ){m90.L2 = m90.L2+1; rej90.L2 = 1}
  if( any(abs(obs.L2 - true.val) > w.95.L2) ){m95.L2 = m95.L2+1; rej95.L2 = 1}
  if( any(abs(obs.L2 - true.val) > w.99.L2) ){m99.L2 = m99.L2+1; rej99.L2 = 1}
  
  
  # fixed L2-norm bandwidth
  rej90.fix = 0; rej95.fix  = 0; rej99.fix  = 0
  if( any(abs(obs.fix - true.val) > w.90.fix) ){m90.fix = m90.fix+1; rej90.fix = 1}
  if( any(abs(obs.fix - true.val) > w.95.fix) ){m95.fix = m95.fix+1; rej95.fix = 1}
  if( any(abs(obs.fix - true.val) > w.99.fix) ){m99.fix = m99.fix+1; rej99.fix = 1}
  
  # log(n^(-0.2)) * fixed L2-norm bandwidth
  rej90.scaleL2 = 0; rej95.scaleL2  = 0; rej99.scaleL2  = 0
  if( any(abs(obs.scaleL2 - true.val) > w.90.scaleL2) ){m90.scaleL2 = m90.scaleL2+1; rej90.scaleL2 = 1}
  if( any(abs(obs.scaleL2 - true.val) > w.95.scaleL2) ){m95.scaleL2 = m95.scaleL2+1; rej95.scaleL2 = 1}
  if( any(abs(obs.scaleL2 - true.val) > w.99.scaleL2) ){m99.scaleL2 = m99.scaleL2+1; rej99.scaleL2 = 1}
  
  
  
  rej90.L2.list = c(rej90.L2.list, rej90.L2); rej90.fix.list = c(rej90.fix.list, rej90.fix); rej90.scaleL2.list = c(rej90.scaleL2.list, rej90.scaleL2)
  
  rej95.L2.list = c(rej95.L2.list, rej95.L2); rej95.fix.list = c(rej95.fix.list, rej95.fix); rej95.scaleL2.list = c(rej95.scaleL2.list, rej95.scaleL2)
  
  rej99.L2.list = c(rej99.L2.list, rej99.L2); rej99.fix.list = c(rej99.fix.list, rej99.fix);  rej99.scaleL2.list = c(rej99.scaleL2.list, rej99.scaleL2)
  
  
  
  if(!delta_constant){
    
    # L2-norm bandwidth
    if( any(abs(obs.L2 - 0) > w.90.L2) ){power90.L2 = power90.L2+1}
    if( any(abs(obs.L2 - 0) > w.95.L2) ){power95.L2 = power95.L2+1}
    if( any(abs(obs.L2 - 0) > w.99.L2) ){power99.L2 = power99.L2+1}
    
    # fixed L2-norm bandwidth
    if( any(abs(obs.fix - 0) > w.90.fix) ){power90.fix = power90.fix+1}
    if( any(abs(obs.fix - 0) > w.95.fix) ){power95.fix = power95.fix+1}
    if( any(abs(obs.fix - 0) > w.99.fix) ){power99.fix = power99.fix+1}
    
    # log(n^(-0.2)) * fixed L2-norm bandwidth
    if( any(abs(obs.scaleL2 - 0) > w.90.scaleL2) ){power90.scaleL2 = power90.scaleL2+1}
    if( any(abs(obs.scaleL2 - 0) > w.95.scaleL2) ){power95.scaleL2 = power95.scaleL2+1}
    if( any(abs(obs.scaleL2 - 0) > w.99.scaleL2) ){power99.scaleL2 = power99.scaleL2+1}
    
   
  }
  
  if(rep %% 10 == 0){print(rep)}
  
  emp.stat.L2       = c(emp.stat.L2,        sqrt(N[j])*max(abs((obs.L2       - true.val)/sqrt(kernel.diag))))
  emp.stat.fix      = c(emp.stat.fix,       sqrt(N[j])*max(abs((obs.fix      - true.val)/sqrt(kernel.diag))))
  emp.stat.scaleL2  = c(emp.stat.scaleL2,   sqrt(N[j])*max(abs((obs.scaleL2  - true.val)/sqrt(kernel.diag))))
}

plan(sequential)




###################
### Saving data ### 
###################

file  = paste0("Simulation/Coverage and power/", folder, "/p ", P[i], "/dep_", N[j], "_rep_", Repitition, "_results.rds")

results = list(

  
  coverage.L2 =  list(q90 = (Repitition - m90.L2)/Repitition,
                      q95 = (Repitition - m95.L2)/Repitition,
                      q99 = (Repitition - m99.L2)/Repitition),
  
  coverage.fix = list(q90 = (Repitition - m90.fix)/Repitition,
                      q95 = (Repitition - m95.fix)/Repitition,
                      q99 = (Repitition - m99.fix)/Repitition),
  
  coverage.scaleL2 = list(q90 = (Repitition - m90.scaleL2)/Repitition,
                          q95 = (Repitition - m95.scaleL2)/Repitition,
                          q99 = (Repitition - m99.scaleL2)/Repitition),
  

  
  power.L2    = if(!delta_constant) list(p90 = power90.L2/Repitition,
                                         p95 = power95.L2/Repitition,
                                         p99 = power99.L2/Repitition) else NULL,
  
  power.fix   = if(!delta_constant) list(p90 = power90.fix/Repitition,
                                         p95 = power95.fix/Repitition,
                                         p99 = power99.fix/Repitition) else NULL,
  
  power.scaleL2   = if(!delta_constant) list(p90 = power90.scaleL2/Repitition,
                                             p95 = power95.scaleL2/Repitition,
                                             p99 = power99.scaleL2/Repitition) else NULL,

  
  q.list.L2        =  list(q90list   = q90.L2.list,     q95list   = q95.L2.list,     q99list   = q99.L2.list),
  q.list.fix        = list(q90list   = q90.fix.list,    q95list   = q95.fix.list,    q99list   = q99.fix.list),
  q.list.scaleL2      = list(q90list = q90.scaleL2.list, q95list  = q95.scaleL2.list, q99list  = q99.scaleL2.list),
  
  sample.list = list( sample.L2        = emp.stat.L2,
                      sample.fix       = emp.stat.fix,
                      sample.scaleL2   = emp.stat.scaleL2),
  
  h.list = data.frame( h.s.L2   = h.s.list.L2,
                       h.d.L2   = h.d.list.L2,
                       
                       rej.L2.90      =  rej90.L2.list,
                       rej.fix.90     =  rej90.fix.list,
                       rej.scaleL2.90 =  rej90.scaleL2.list,
                       
                       rej.L2.95      =  rej95.L2.list,  
                       rej.fix.95     =  rej95.fix.list, 
                       rej.scaleL2.95 =  rej95.scaleL2.list,
                       
                       rej.L2.99      =  rej99.L2.list,
                       rej.fix.99     =  rej99.fix.list,
                       rej.scaleL2.99 =  rej99.scaleL2.list
                       
                 )
)

saveRDS(results,  file)









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


load("Simulation/Compare_independent_dependent_Multiplier_Bootstrap.RData")

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

# Setting Parameters 
C_inv  = 1.2
n      = 400
nd     = n * C_inv
p      = 25
pd     = 100
p.eval = 75


#Gernate data sets
if(delta_constant){set.seed(1144)}else{set.seed(1169)}
x.s.design = (1:p - 0.5)/p
Y.s        = t(mu(x.s.design) + 
                 t(sim.d.OU(n, t = x.s.design,rho_B = 0.5,tau = 0,sigma = 4)) + 
                 matrix(rnorm(length(x.s.design)*n, 0, 0.1),length(x.s.design), n))

if(delta_constant){set.seed(1148)}else{set.seed(1179)}
xd.design = (1:pd - 0.5)/pd
Y.d       = t(mu_d(xd.design,delta_constant) + 
                t(sim.d.OU(nd, t = xd.design,rho_B = 0.5,tau = 0,sigma = 4)) + 
                matrix(rnorm(length(xd.design)*nd, 0, 0.1),length(xd.design), nd))

#########################################
# Cumputing long run covariance kernel  #
#########################################

#sparse 
w       = readRDS(paste0("Simulation/weights/p ",p,"/full_w_lag0_",n,".rds"))
Gamma.s = eval.weights(w,observation.transformation(Y.s,grid.type = "less"))


# Optional: For parallelizing computations ###############
options(future.globals.maxSize = 20 * 1024^3)           ##  
plan(multisession, workers = future::availableCores()-2)##

Gamma.s.part2 = Reduce(`+`,future_lapply(1:3, function(k) {
  
  w.lag = readRDS(paste0("Simulation/weights/p ",p,"/full_w_lag",k,"_",n,".rds"))
  lag.Gamma.s = eval.weights(w.lag, observation.transformation(Y.s, lag = k, grid.type = "lesseq"), lag = k)
  lag.Gamma.s = (1 - k/(3 + 1) ) * (lag.Gamma.s + t(lag.Gamma.s))
  
  lag.Gamma.s 
}))

plan(sequential)


Gamma.s   = Gamma.s + Gamma.s.part2
P.Gamma.s = diag(Gamma.s) - rowMeans(Gamma.s) - colMeans(Gamma.s) + mean(Gamma.s) 


#dense 
wd      = readRDS(paste0("Simulation/weights/pd 100 ",delta_constant,"/full_w_lag0_",n,".rds"))
Gamma.d = eval.weights(wd,observation.transformation(Y.d,grid.type = "less"))


# Optional: For parallelizing computations ###############
options(future.globals.maxSize = 20 * 1024^3)           ##  
plan(multisession, workers = future::availableCores()-2)##

Gamma.d.part2 = Reduce(`+`,future_lapply(1:3, function(k) {
  
  wd.lag = readRDS(paste0("Simulation/weights/pd 100 ",delta_constant,"/full_w_lag",k,"_",n,".rds"))
  lag.Gamma.d = eval.weights(wd.lag, observation.transformation(Y.d, lag = k, grid.type = "lesseq"), lag = k)
  lag.Gamma.d = (1 - k/(3 + 1) ) * (lag.Gamma.d + t(lag.Gamma.d))
  
  lag.Gamma.d 
}))

plan(sequential)


Gamma.d   = Gamma.d + Gamma.d.part2
P.Gamma.d = diag(Gamma.d) - rowMeans(Gamma.d) - colMeans(Gamma.d) + mean(Gamma.d) 


### long run kernel ################################# 
Gamma   = diag(Gamma.s) + n/nd* diag(Gamma.d)       #
P.Gamma = P.Gamma.s     + n/nd* P.Gamma.d           # 
#####################################################

rm(w,wd,Gamma.s,Gamma.d, Gamma.s.part2, Gamma.d.part2, P.Gamma.s, P.Gamma.d)

# Mean and Difference Estimation #
bandwidth  = readRDS(paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds"))
estimation = est.results(Y.s,Y.d,bandwidth[c(2,4),11], p.eval = p.eval)



#################################
# Applying Multiplier Bootstrap #
#################################

# Optional: For parallelizing computations ###############
options(future.globals.maxSize = 20 * 1024^3)           ##  
plan(multisession, workers = future::availableCores()-2)##

dep.MB.result = q.MB(Y.s, Y.d, estimation$sparse, estimation$dense, P.Gamma, (1:p.eval)/p.eval, bandwidth[c(2,4),11], alpha = 0.95, B = 1000, depend = T, int = T, H0 = 0)
ind.MB.result = q.MB(Y.s, Y.d, estimation$sparse, estimation$dense, P.Gamma, (1:p.eval)/p.eval, bandwidth[c(2,4),11], alpha = 0.95, B = 1000, depend = F, int = T, H0 = 0)

plan(sequential)


delta.centered     = delta((1:p.eval-0.5)/p.eval, constant = delta_constant)-integrate(function(x){delta(x,constant = delta_constant)},0,1)$value 
delta.centered.est = estimation$delta$ESTIMATE - estimation$delta_int$ESTIMATE

delta.centered.df = data.frame(x = (1:p.eval-0.5)/p.eval, true.val = delta.centered, est.val = delta.centered.est)

KI.dep = data.frame(x = (1:p.eval-0.5)/p.eval, UP = delta.centered.est + sqrt(P.Gamma/(n))* dep.MB.result$quantile,LOW = delta.centered.est - sqrt(P.Gamma/(n))* dep.MB.result$quantile)
KI.ind = data.frame(x = (1:p.eval-0.5)/p.eval, UP = delta.centered.est + sqrt(P.Gamma/(n))* ind.MB.result$quantile,LOW = delta.centered.est - sqrt(P.Gamma/(n))* ind.MB.result$quantile)

ggplot() +
  geom_ribbon(aes(x = x, ymin = LOW, ymax = UP, fill = "dependent"),data = KI.dep, alpha = 0.2, col = NA) +
  geom_ribbon(aes(x = x, ymin = LOW, ymax = UP, fill = "independent"),data = KI.ind, col = "black", alpha = 0,size = 0.7, lty = 2)+
  geom_line(aes(x, est.val, colour = "hat Difference"), data = delta.centered.df, lty = 1, size = 0.8) +
  geom_line(aes(x, true.val, colour = "Difference"), data = delta.centered.df, lty = 1, size = 0.6) +
  scale_fill_manual("Confidence bands:", values = c("dependent" = "black", "independent" = "#F8766D"),
                    labels = c("dependent" = expression("with "*hat(q)[0.95]^{DMB}),
                               "independent" = expression("with "*hat(q)[0.95]^{IMB})),
                    breaks = c("dependent","independent")) +
  
  scale_colour_manual("Curves:", values = c("hat Difference" = "yellow","Difference" = "#F8766D"),
                      labels = c("hat Difference" = expression(hat(delta) - integral(hat(delta), "")*" d"*lambda),
                                 "Difference" = expression(delta - integral(delta, "")*" d"*lambda))) +
  theme(plot.subtitle = element_text(size =14),
        legend.text = element_text(size =12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 15),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 13),     
        axis.text.y  = element_text(size = 13))+
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(0.3, "cm"))+
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    colour = guide_legend(nrow = 2, byrow = TRUE)
  )

ggsave("Simulation/pictures/Comparison_CB.png", width = 26, height = 14, units = "cm", dpi = 300)

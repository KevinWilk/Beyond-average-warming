library(MASS)
library(ggplot2)
library(reshape2)
library(locpol)
library(interp)
library(stats)
library(future.apply)
library(parallel)
library(cowplot)
library(biLocPol)
library(plotly)
library(tidyr)
library(goffda)
library(crayon)
library(future)
library(tidyverse)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)


source("Mod_biLocPol.R")

source("Simulation/Simulation_Mod_biLocPol.R") #Adapted to grid from Berger and Holzmann (2025)
source("Simulation/function_data_generating.R")
source("Simulation/function_bandwidth_selection.R")
source("Simulation/function_evaluation.R")









###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
#########                                                                                ##################################
######### Simulation illustration for n = 25, p = 25 and \tilde{n} = 30, \tilde{p} = 100 ##################################
#########                                                                                ##################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################





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



# Setting Parameters 
x = seq(0, 1, 0.01)
n  = 25
nd = n*1.2
p  = 25
pd = 100
obs.s = (1:p  - 0.5)/p
obs.d = (1:pd - 0.5)/pd


if(delta_constant){set.seed(1144)}else{set.seed(1167)}
#### target sample simulation #####################################################
Y.s = mu(obs.s) + #################################################################
      t(sim.d.OU(n, t = obs.s,rho_B = 0.5,tau = 0,sigma = 4)) + ###################
      matrix(rnorm(length(obs.s)*n, 0, 0.1),length(obs.s), n)######################
mu.s.df = data.frame(x = x, y = mu(x), fac = "mu")#################################
Y.s.df = melt(data.frame(obs.s, Y.s), id.vars = "obs.s", ##########################
              variable.name = "i", value.name = "Y")###############################
Y.s.mean = data.frame(obs.s, Y.s= rowMeans(Y.s))###################################
###################################################################################

if(delta_constant){set.seed(1148)}else{set.seed(1177)}
#### source sample simulation ###################################################################
Y.d = mu_d(obs.d,delta_constant) + ##############################################################
                t(sim.d.OU(nd, t = obs.d,rho_B = 0.5,tau = 0,sigma = 4)) + ######################
                matrix(rnorm(length(obs.d)*nd, 0, 0.1),length(obs.d), nd)########################
#################################################################################################
mu.d.df = data.frame(x = x, Y.d = mu_d(x,delta_constant), fac = "mu_d")##########################
Y.d.df = melt(data.frame(obs.d, Y.d), id.vars = "obs.d", ########################################
              variable.name = "i", value.name = "Y.d")###########################################
Y.d.mean = data.frame(obs.d, Y.d = rowMeans(Y.d))################################################
#################################################################################################


#### Bandwidth for source sample and target ##########################################################################
bandwidth  = readRDS(paste0("Simulation/optimal bandwidth/Results/delta_constant_",delta_constant,"_Bandwidth.rds"))##
bw.d  = bandwidth[2,2]; bw.delta  = bandwidth[4,2] ###################################################################
######################################################################################################################



#### Estimation on grid x #####################################################
Y.d.est.x = locPolSmootherC(obs.d, Y.d.mean$Y.d, x, bw.d, 2, EpaK)#############
Y.d.est.x = data.frame(x = x, Y.d.hat = Y.d.est.x$beta0)#######################
###############################################################################


#### Estimation on grid obs.s ################################################
Y.d.est = locPolSmootherC(obs.d, Y.d.mean$Y.d, obs.s, bw.d, 2, EpaK)##########
Y.d.est.obs.s = data.frame(x = obs.s, Y.d.hat = Y.d.est$beta0)################
##############################################################################


#### Creating Residuals and Means on grid obs.s ######################
res = Y.s-Y.d.est.obs.s$Y.d.hat#######################################
res.mean = data.frame(obs.s, Y = rowMeans(res))#######################
######################################################################


#### Estimation mu grid x  with TF #########################################################
delta.est = locPolSmootherC(obs.s, Y.s.mean$Y-Y.d.est.obs.s$Y.d.hat, x, bw.delta, 2, EpaK)##
mu.est    = locPolSmootherC(obs.s, Y.s.mean$Y, x, bw.delta, 2, EpaK) #######################
Y.d.est = locPolSmootherC(obs.d, Y.d.mean$Y.d, x, bw.d, 2, EpaK)############################
Y.res.est = data.frame(x = x, Y.res.hat = delta.est$beta0+Y.d.est$beta0)####################
############################################################################################


ggplot() + 
  geom_point(aes(obs.s, Y, col = i), Y.s.df, size = 2,alpha = .5)+
  geom_line(aes(obs.s, Y, col = i), Y.s.df, alpha = .5, size = 1, lty = 2, show.legend = F)+
  geom_line(aes(x, y), mu.s.df,size = 1.5) +
  geom_line(aes(x, Y.res.hat), Y.res.est, lty = 2,col = "red",size = 2) +
  theme(plot.subtitle = element_text(size =26),
        legend.position = "none",
        legend.text = element_text(size =17),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 17),     
        axis.title.y = element_text(size = 17),
        axis.text.x  = element_text(size = 24),     
        axis.text.y  = element_text(size = 24))+
  ylim(-12.5,11)+
  labs(y=NULL, x = NULL, 
       subtitle = bquote("Simulation setup of " * mu^{"[s]"} * ": n = " * .(n) * ", p = " * .(p))) 

ggsave(paste0("Simulation/Figures/constant ",delta_constant," simulation setup sparse.png"), width = 18, height = 17, units = "cm", dpi = 300)


ggplot() + 
  geom_line(aes(obs.d, Y.d, col = i), Y.d.df, alpha = .3, size = 1, lty = 2, show.legend = F)+
  geom_point(aes(obs.d, Y.d, col = i), Y.d.df, size = 1.5,alpha = .3)+
  geom_line(aes(x, Y.d), mu.d.df,size = 1.5) +
  geom_line(aes(x, Y.d.hat), Y.d.est.obs.s, lty = 2,col = "red",size = 2) +
  theme(plot.subtitle = element_text(size =26),
        legend.position = "none",
        legend.text = element_text(size  = 17),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 17),     
        axis.title.y = element_text(size = 17),
        axis.text.x  = element_text(size = 24),     
        axis.text.y  = element_text(size = 24))+
  ylim(-12.5,11)+
  labs(y=NULL, x = NULL, 
       subtitle = bquote("Simulation setup of " * mu^{"[d]"} * ": "*tilde(n)*" = " * .(nd) * ", "*tilde(p)*" = " * .(pd))) 

ggsave(paste0("Simulation/Figures/constant ",delta_constant," simulation setup dense.png"), width = 18, height = 17, units = "cm", dpi = 300)



delta.est.df = data.frame(x = x, delta.hat = delta.est$beta0)
delta.sample.df   = melt(data.frame(obs.s, res), id.vars = "obs.s", variable.name = "i", value.name = "Y_Y.d")
delta_data = data.frame(x = x, y = delta(x,delta_constant), fac = "delta")

ggplot() + 
  geom_line(aes(obs.s, Y_Y.d, col = i), delta.sample.df, alpha = .5, size = 1, lty = 2, show.legend = F)+
  geom_point(aes(obs.s, Y_Y.d, col = i), delta.sample.df, size = 2,alpha = .4)+
  geom_line(aes(x, y), delta_data,size = 1.5) +
  geom_line(aes(x, delta.hat), delta.est.df, lty = 2,col = "red",size = 2) +
  theme(plot.subtitle = element_text(size =26),
        legend.position = "none",
        legend.text = element_text(size =17),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 17),     
        axis.title.y = element_text(size = 17),
        axis.text.x  = element_text(size = 24),     
        axis.text.y  = element_text(size = 24))+
  ylim(-12,11)+
  labs(y=NULL, x = NULL, 
       subtitle = bquote("Simulation setup of " * delta * " = " * mu^{"[s]"}-  mu^{"[d]"} )) 

ggsave(paste0("Simulation/Figures/constant ",delta_constant," simulation setup delta.png"), width = 18, height = 17, units = "cm", dpi = 300)




















#####################################################################################################
#####################################################################################################
#####################################################################################################
#########                                                          ##################################
#########   Boxplots of hv cross-validation in a 5-fold framework  ##################################
#########                                                          ##################################
#####################################################################################################
#####################################################################################################
#####################################################################################################






# Setting Parameters 
P  = c(25,50,75)
pd = 100


lag_Gamma_s_25_sample  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_25_sample.rds")
lag_Gamma_s_50_sample  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_50_sample.rds")
lag_Gamma_s_75_sample  = readRDS("Simulation/optimal bandwidth/Results/lag_Gamma_s_75_sample.rds")
lag_Gamma_d_100_sample = readRDS(paste0("Simulation/optimal bandwidth/Results/",if(delta_constant){"Delta = 0"}else{"Delta != 0"},"/lag_Gamma_d_100_sample.rds"))



bw  = c(c(t(lag_Gamma_s_25_sample[[11]])),c(t(lag_Gamma_s_50_sample[[11]])),c(t(lag_Gamma_s_75_sample[[11]])),c(t(lag_Gamma_d_100_sample[[11]])))

lag = rep(rep(0:4, each = 1000), times = 4)


results.hv = data.frame(bandwidth = bw, lag = factor(lag), p = factor(rep(rep(c(P,pd)), each = 5000)))

ggplot(results.hv |> filter(lag %in% 0:3), aes(x = lag, y = bandwidth, color = p)) +
  geom_boxplot(position = position_dodge(width = 0.8), size = 1) +
  labs(x = ~Lag, y = "Bandwidth") +
  lims(y = c(0,1)) +
  theme(plot.title = element_text(size =26),
        legend.text = element_text(size =17),
        strip.text = element_text(size = 0),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),     
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 17),     
        axis.text.y  = element_text(size = 17))+
  facet_wrap(~ p, ncol = 5)


ggsave("Simulation/Figures/hv_bandwidth_selection.png", width = 26, height = 14, units = "cm", dpi = 300)

















#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#########                                                                                  ##################################
#########  Coparison of dependent with independent Multiplier Bootstrap under Alternative  ##################################
#########                                                                                  ##################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################




load("Simulation/Figures/Compare_independent_dependent_Multiplier_Bootstrap.RData")

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


delta.centered     = delta((1:p.eval-0.5)/p.eval, 
                           constant = delta_constant)-integrate(function(x){delta(x,constant = delta_constant)},0,1)$value 
delta.centered.est = estimation$delta$ESTIMATE - estimation$delta_int$ESTIMATE

delta.centered.df = data.frame(x = (1:p.eval-0.5)/p.eval, true.val = delta.centered, est.val = delta.centered.est)

KI.dep = data.frame(x = (1:p.eval-0.5)/p.eval, 
                    UP = delta.centered.est + sqrt(P.Gamma/(n))* dep.MB.result$quantile,
                    LOW = delta.centered.est - sqrt(P.Gamma/(n))* dep.MB.result$quantile)

KI.ind = data.frame(x = (1:p.eval-0.5)/p.eval, 
                    UP = delta.centered.est + sqrt(P.Gamma/(n))* ind.MB.result$quantile,
                    LOW = delta.centered.est - sqrt(P.Gamma/(n))* ind.MB.result$quantile)

ggplot() +
  geom_ribbon(aes(x = x, ymin = LOW, ymax = UP, fill = "dependent"),data = KI.dep, alpha = 0.2, col = NA) +
  geom_ribbon(aes(x = x, ymin = LOW, ymax = UP, fill = "independent"),data = KI.ind, col = "black", alpha = 0,size = 1.2, lty = 2)+
  geom_line(aes(x, est.val, colour = "hat Difference"), data = delta.centered.df, lty = 2, size = 1.2) +
  geom_line(aes(x, true.val, colour = "Difference"), data = delta.centered.df, lty = 1, size = 1) +
  
  labs(subtitle = bquote("n = " * .(n) * ", p = " * .(p) *", " *tilde(n)*" = " * .(nd) * " and "*tilde(p)*" = " * .(pd)),
       x = NULL,
       y = NULL) +
  
  scale_fill_manual("Confidence bands:", 
                    values = c("dependent" = "black", "independent" = "black"),
                    labels = c("dependent" = expression("with "*hat(q)[0.95]^{DMB}),
                               "independent" = expression("with "*hat(q)[0.95]^{IMB})),
                    breaks = c("dependent","independent")) +
  
  scale_colour_manual("Curves:", 
                      values = c("hat Difference" = "red","Difference" = "black"),
                      labels = c("hat Difference" = expression(hat(delta) - integral(hat(delta), "")*" d"*lambda),
                                 "Difference" = expression(delta - integral(delta, "")*" d"*lambda))) +
  
  theme(plot.subtitle = element_text(size = 24),
        legend.text = element_text(size  = 22),
        legend.title = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 22),     
        axis.text.y  = element_text(size = 22),
        
        legend.position   = "right",
        legend.box        = "vertical",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width  = unit(1.1, "cm")) +
  
  guides(fill   = guide_legend(nrow = 2, byrow = TRUE),
         colour = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Simulation/Figures/Comparison_CB.png", width = 30, height = 14, units = "cm", dpi = 300)




















##########################################################################
##########################################################################
####                                     #################################
####  Illustation of coverage and power  #################################
####                                     #################################
##########################################################################
##########################################################################


# Parameters
N   = c(25,50,100,150,200,250,300,350,400)
P   = c(25,50,75)
Rep = 1000


delta_constant = FALSE
integral       = FALSE

##########################################
if(delta_constant){                  #####
  if(integral){                      #####
    folder = "Under H0"              #####
  }else{folder = "H0 not centered"}  #####
  ########################################
}else{                               #####
  if(integral){                      #####
    folder = "Under H1"              #####
  }else{folder = "H1 not centered"}  #####
}                                    #####
##########################################



###################
# Under H0 or H1 ##
###################

dep.value.c = c()
ind.value.c = c()

for(i in 1:length(P)){
  
  dep.coverage.90  = c()
  dep.coverage.95  = c()
  dep.coverage.99  = c()
  
  ind.coverage.90  = c()
  ind.coverage.95  = c()
  ind.coverage.99  = c()
  
  for(j in 1:9){
    
    dep.c90 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_coverage90.rds"))
    dep.c95 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_coverage95.rds"))
    dep.c99 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_coverage99.rds"))
    
    dep.coverage.90  = c(dep.coverage.90,dep.c90)
    dep.coverage.95  = c(dep.coverage.95,dep.c95)
    dep.coverage.99  = c(dep.coverage.99,dep.c99)
    
    
    #ind.c90 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_coverage90.rds"))
    #ind.c95 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_coverage95.rds"))
    #ind.c99 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_coverage99.rds"))
    
    #ind.coverage.90  = c(ind.coverage.90,ind.c90)
    #ind.coverage.95  = c(ind.coverage.95,ind.c95)
    #ind.coverage.99  = c(ind.coverage.99,ind.c99)
    
  }
  
  dep.value.c = c(dep.value.c,dep.coverage.90,dep.coverage.95,dep.coverage.99)
  
  #ind.value.c = c(ind.value.c,ind.coverage.90,ind.coverage.95,ind.coverage.99)
  
}

dep.coverage.df    = data.frame(p = rep(paste0("p = ", P), each = length(N)*3) ,n = rep(N, times = length(P)*3), quantile.est = rep(c("90%","95%","99%"), each = length(N)), emp.coverage = dep.value.c, type = "i = DMB")
dep.endpoints      = dep.coverage.df |> group_by(p,quantile.est) |> filter(emp.coverage == max(emp.coverage))
ind.coverage.df    = data.frame(p = rep(paste0("p = ", P), each = length(N)*3) ,n = rep(N, times = length(P)*3), quantile.est = rep(c("90%","95%","99%"), each = length(N)), emp.coverage = ind.value.c, type = "i = IMB")
coverage.df        = rbind(dep.coverage.df,ind.coverage.df)

##########################################
if(delta_constant){                  #####
  if(integral){                      #####
    folder = "H0 centered"           #####
  }else{folder = "H0 not centered"}  #####
  ########################################
}else{                               #####
  if(integral){                      #####
    folder = "H1 centered"           #####
  }else{folder = "H1 not centered"}  #####
}                                    #####
##########################################


ggplot() + 
  geom_line( aes(x = n, y = emp.coverage, colour = quantile.est, linetype = type), data = dep.coverage.df, size = 1.6) +
  geom_text(data = dep.endpoints, aes(x = 400, y = emp.coverage, label = paste0(round(emp.coverage*100, digits = 1),"%"), colour = quantile.est), hjust = -0.15, size = 6, show.legend = FALSE) +
  
  scale_colour_manual(name = expression(1-alpha),values = c("90%" = "#F8766D", "95%" = "#7CAE00", "99%" = "#00BFC4") ) +
  
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text = element_text(size = 25),
        plot.margin = margin(5.5, 27, 5.5, 5.5),     
        legend.box.margin = margin(0, 3, 0, 0),     
        legend.title = element_text(size = 25),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        
        legend.key.width  = unit(1.1, "cm")) +
  labs(linetype = expression("with "*{hat(q)^{i}}[1-alpha])) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0,1)) +
  xlim(25,475) +
  labs(y = "empirical coverage", x = "n", title = bquote(tilde(n)* " = 1.2 n  and  " *tilde(p)* " = 100")) +
  facet_grid(~factor(p))

ggsave(paste0("Simulation/Figures/coverage rate ",folder,".png"), width = 39, height = 16, units = "cm", dpi = 300)










####################
# Under H1: power ##
####################

# Parameters
N   = c(25,50,75,100,150,200,250,300,350,400)

dep.value.p = c()
ind.value.p = c()

for(i in 1:length(P)){
  
  dep.power.90  = c()
  dep.power.95  = c()
  dep.power.99  = c()
  
  ind.power.90  = c()
  ind.power.95  = c()
  ind.power.99  = c()
  
  for(j in 1:10){
    
    dep.p90 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_power90.rds"))
    dep.p95 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_power95.rds"))
    dep.p99 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/dep_",N[j],"_rep_",Rep,"_power99.rds"))
    
    dep.power.90  = c(dep.power.90,dep.p90)
    dep.power.95  = c(dep.power.95,dep.p95)
    dep.power.99  = c(dep.power.99,dep.p99)
    
    
    #ind.p90 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_power90.rds"))
    #ind.p95 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_power95.rds"))
    #ind.p99 = readRDS(paste0("Simulation/Coverage and power/",folder,"/p ",P[i],"/ind_",N[j],"_rep_",Rep,"_power99.rds"))
    
    #ind.power.90  = c(ind.power.90,ind.p90)
    #ind.power.95  = c(ind.power.95,ind.p95)
    #ind.power.99  = c(ind.power.99,ind.p99)
    
  }
  
  dep.value.p = c(dep.value.p,dep.power.90,dep.power.95,dep.power.99)
  
  #ind.value.p = c(ind.value.p,ind.power.90,ind.power.95,ind.power.99)
  
}

dep.power.df    = data.frame(p = rep(paste0("p = ", P), each = length(N)*3) ,n = rep(N, times = length(P)*3), quantile.est = rep(c("90%","95%","99%"), each = length(N)), emp.power = dep.value.p, type = "i = DMB")


##########################################
if(delta_constant){                  #####
  if(integral){                      #####
    folder = "H0 centered"           #####
  }else{folder = "H0 not centered"}  #####
  ########################################
}else{                               #####
  if(integral){                      #####
    folder = "H1 centered"           #####
  }else{folder = "H1 not centered"}  #####
}                                    #####
##########################################

ggplot() + 
  geom_line( aes(x = n, y = emp.power, colour = quantile.est, linetype = type), data = dep.power.df, size = 1.6) +

  scale_colour_manual(name = expression(1-alpha),values = c("90%" = "#F8766D", "95%" = "#7CAE00", "99%" = "#00BFC4") ) +
  
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text = element_text(size = 25),
        plot.margin = margin(5.5, 27, 5.5, 5.5),     
        legend.box.margin = margin(0, 3, 0, 0),     
        legend.title = element_text(size = 25),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        
        legend.key.width  = unit(1.1, "cm")) +
  labs(linetype = expression("with "*{hat(q)^{i}}[1-alpha])) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0,1)) +
  xlim(25,475) +
  labs(y = "power", x = "n", title = "") +
  facet_grid(~factor(p))


ggsave(paste0("Simulation/Figures/power ",folder,".png"), width = 39, height = 16, units = "cm", dpi = 300)










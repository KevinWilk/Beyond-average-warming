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
            variable.name = "i", value.name = "Y")#################################
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
             variable.name = "i", value.name = "Y.d")############################################
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

ggsave(paste0("Simulation/pictures/constant ",delta_constant," simulation setup sparse.png"), width = 18, height = 17, units = "cm", dpi = 300)


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

ggsave(paste0("Simulation/pictures/constant ",delta_constant," simulation setup dense.png"), width = 18, height = 17, units = "cm", dpi = 300)



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

ggsave(paste0("Simulation/pictures/constant ",delta_constant," simulation setup delta.png"), width = 18, height = 17, units = "cm", dpi = 300)




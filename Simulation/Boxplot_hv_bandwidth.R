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
delta_constant = TRUE     ################
                          ################
##########################################



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
















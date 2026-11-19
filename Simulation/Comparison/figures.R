library(crayon)
library(ggplot2)
library(reshape2)
library(locpol)
library(biLocPol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(lubridate)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)


source("Mod_biLocPol.R")
source("weather temperature/functions.R")
source("weather temperature/bandwidth selection/function_bandwidth_selection.R")


# Parameters
N   = c(25,50,75,100,150,200,250,300,350,400,600,1000)
P   = c(25,50,75,15,16,18,20)
Rep = 1000


delta_constant = FALSE

##########################################
if(delta_constant){                  #####
  folder = "H0 (centered)"           #####
  name   = "H0"                      #####
    ######################################
}else{                               #####
  folder = "H1 (centered)"           #####
  name   = "H1"                      #####
}                                    #####
##########################################






value.reject = c()

for(i in c(4:7,1:3)){
  
  reject.90  = c()
  reject.95  = c()
  reject.99  = c()
  
  
  for(j in c(1:12)){
    
    reject = readRDS(paste0("Simulation/Comparison/",folder,"/p ",P[i],"/dep_",N[j],"_rep_1000_reject_list.rds"))
    
    reject90 = reject$rate90
    reject95 = reject$rate95
    reject99 = reject$rate99
    
    reject.90  = c(reject.90,reject90)
    reject.95  = c(reject.95,reject95)
    reject.99  = c(reject.99,reject99)
    
  }
  
  value.reject = c(value.reject,reject.90,reject.95,reject.99)
  
}

reject.df = data.frame(p = rep(paste0("p = ", P[c(4:7,1:3)]), each = length(N)*3) ,
                       n = rep(N, times = length(P[c(4:7,1:3)]*3)), 
                       quantile.est = rep(c("90%","95%","99%"), each = length(N)), 
                       emp.rate = value.reject)


endpoints99      = reject.df |> group_by(p,quantile.est) |> slice_max(n, n = 1, with_ties = FALSE) |> filter(quantile.est == "99%") |> filter( !(p %in% c("p = 25", "p = 50", "p = 75")) )
endpoints95      = reject.df |> group_by(p,quantile.est) |> slice_max(n, n = 1, with_ties = FALSE) |> filter(quantile.est == "95%") |> filter( !(p %in% c("p = 25", "p = 50", "p = 75")) )
endpoints90      = reject.df |> group_by(p,quantile.est) |> slice_max(n, n = 1, with_ties = FALSE) |> filter(quantile.est == "90%")
endpoints = rbind(endpoints99,endpoints95,endpoints90)



ggplot() + 
    
  geom_line( aes(x = n, y = emp.rate, colour = quantile.est, linetype = quantile.est), data = reject.df, size = 1.6) +
  geom_text(data = endpoints, aes(x = 1000, y = emp.rate, label = paste0(round(emp.rate*100, digits = 1),"%"), colour = quantile.est), hjust = -0.15, size = 9, show.legend = FALSE) +
  
  scale_colour_manual(name = expression("with "*{hat(q)^{DMB}}[1-alpha]*":"),
                      breaks = c("90%", "95%", "99%"),
                      values = c("90%" = "#F8766D","95%" = "#7CAE00","99%" = "#00BFC4")) +
  
  scale_linetype_manual(values = c("90%" = 1,"95%" = 2,"99%" = 4), guide = "none") +
  
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 4),linewidth = 1.6))) +
  
  theme(plot.title = element_text(size = 28),
        plot.subtitle = element_text(size = 28),
        legend.text = element_text(size = 28),
        strip.text = element_text(size = 28),
        plot.margin = margin(5.5, 27, 5.5, 5.5),     
        legend.box.margin = margin(0, 3, 0, 0),     
        legend.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),     
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 28),     
        axis.text.y  = element_text(size = 28),
        
        legend.key.width  = unit(1.1, "cm")) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1)) +
  scale_x_continuous(breaks = c(0, 500, 1000), limits = c(0,1380)) +
  labs(y = "rejection rate", x = "n", title = bquote(tilde(n)* " = 1.2 n  and  " *tilde(p)* " = 100")) +
  facet_grid(. ~ p)

ggsave(paste0("Simulation/Comparison/rejection_rate_ ",name,".png"), width = 82, height = 18, units = "cm", dpi = 300)


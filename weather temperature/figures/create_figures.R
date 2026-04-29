####################################################################################################
####################################################################################################
###############              #######################################################################
###############    Figures   #######################################################################
###############              #######################################################################
####################################################################################################
####################################################################################################


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
library(tidyr)
library(tibble)
library(gridExtra)


source("weather temperature/functions.R")


###############################################################################
k = 1  # Change to: 1 (Berlin)                                               ##
#                   2 (Frankfurt am Main)                                    ##
#                   3 (Hamburg)                                              ##
#                   4 (Munich)                                               ##
data.example     = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")     ##
data.example.pic = list("Berlin", "Frankfurt am Main", "Hamburg", "Munich")  ##
load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))    ##
###############################################################################









#############################################################################################
##### Plot of inference autocovariance test results of each month ###########################
#############################################################################################


# Loading estimated lag covariance kernels for lag = 0,...,12 
cov.s.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_s_",data.example[[k]],".rds")) 
cov.d.list = readRDS(paste0("weather temperature/long run kernel/Results/list_Gamma_d_",data.example[[k]],".rds")) 


start.24 = 21; end.24 = 117 

hat.rho = function(h, cov.data){lapply(1:12, function(m) {cov.lag    = sqrt(mean((cov.data[[m]][[h]][start.24:end.24,start.24:end.24])^2))
cov        = mean(diag( cov.data[[m]][[1]][start.24:end.24,start.24:end.24]))
cov.lag/cov })}

rho.hat.s.df = as.data.frame(do.call(rbind, lapply(2:13, function(k) {unlist(hat.rho(k, cov.s.list))})))
colnames(rho.hat.s.df) = month.name
rownames(rho.hat.s.df) = paste0("lag ", 1:12)
rho.hat.s.df = rho.hat.s.df |> rownames_to_column("lag") |> mutate(lag = as.integer(gsub("lag\\s*", "", lag))) |> pivot_longer(cols = -lag,names_to = "MONTH",values_to = "value") |> mutate(MONTH = factor(MONTH, levels = month.name))

rho.hat.d.df = as.data.frame(do.call(rbind, lapply(2:13, function(k) {unlist(hat.rho(k, cov.d.list))})))
colnames(rho.hat.d.df) = month.name
rownames(rho.hat.d.df) = paste0("lag ", 1:12)
rho.hat.d.df = rho.hat.d.df |> rownames_to_column("lag") |> mutate(lag = as.integer(gsub("lag\\s*", "", lag))) |> pivot_longer(cols = -lag,names_to = "MONTH",values_to = "value") |> mutate(MONTH = factor(MONTH, levels = month.name))



# Loading test.s.val and test.d.val                           
test.s.val = readRDS(paste0("weather temperature/long run kernel/Results/test_max_lag_s_",data.example[[k]],".rds")) 
test.d.val = readRDS(paste0("weather temperature/long run kernel/Results/test_max_lag_d_",data.example[[k]],".rds")) 



ggplot() +
  geom_bar( aes(x = lag,     y = value), stat = "identity", fill = "turquoise",width = 0.4, data = rho.hat.d.df, alpha = 0.7) +
  geom_line(aes(x = lag,     y = value), test.d.val, col = "green3", size = 1.1, lty = 2)+
  geom_line(aes(x = lag+0.3, y = value), test.s.val, col = "blue", size = 1.1, lty = 2)+
  geom_bar( aes(x = lag+0.3, y = value),  stat = "identity", fill = "darkblue", width = 0.4, data = rho.hat.s.df, alpha = 0.5) +
  labs(x = "Lag", y = bquote(hat(rho)[Lag]) ,title = "") +
  
  theme(plot.title   = element_text(size = 15),
        legend.text  = element_text(size = 15),
        strip.text   = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 15),     
        axis.text.y  = element_text(size = 15)) +
  
  scale_y_continuous(breaks = c(0,0.25,0.5), limits = c(0,0.7)) +
  scale_x_continuous(breaks = c(1,4,8,12), limits = c(0.5,12.8)) +
  facet_wrap(~ MONTH)

ggsave(paste0("weather temperature/figures/lag_test_",data.example[[k]],".png"), width = 32, height = 18, units = "cm", dpi = 300)








#############################################################################################
##### Plot of standard deviation (lag = 0) of each month ####################################
#############################################################################################


# Loading maximum lag of long run kernel 
max.lag = readRDS(paste0("weather temperature/long run kernel/Results/max_lag_",data.example[[k]],".rds"))


start.24 = 21; end.24 = 117 

sqrt.d.month = lapply(1:12, function(m) {sqrt(diag(cov.d.list[[m]][[1]][start.24:end.24,start.24:end.24]))})
sqrt.s.month = lapply(1:12, function(m) {sqrt(diag(cov.s.list[[m]][[1]][start.24:end.24,start.24:end.24]))})
lr.sqrt.d.month = lapply(1:12, function(m) {sqrt(Reduce(`+`,lapply(1:(max.lag[m,2]+1),  function(j){diag(cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]) * ifelse(j >= 2, 2*(1-(j-1)/(max.lag[m,2]+1)), 1) })))})
lr.sqrt.s.month = lapply(1:12, function(m) {sqrt(Reduce(`+`,lapply(1:(max.lag[m,2]+1),  function(j){diag(cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]) * ifelse(j >= 2, 2*(1-(j-1)/(max.lag[m,2]+1)), 1) })))})
eval.tibble = tibble(TIME = rep(hms::as_hms(c(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min"), as.POSIXct("1970-01-01 23:59:59"))), times = 12))
sqrt.month = data.frame(TIME = eval.tibble, SD.sparse = unlist(sqrt.s.month), SD.sparse.lr = unlist(lr.sqrt.s.month), SD.dense = unlist(sqrt.d.month) , SD.dense.lr = unlist(lr.sqrt.d.month), MONTH = rep(month.abb,each = 97))
sqrt.month$MONTH = factor(sqrt.month$MONTH,levels = month.abb)
rm(lr.sqrt.s.month,lr.sqrt.d.month,sqrt.d.month,sqrt.s.month,eval.tibble)

month.colors = c("Dec" = "#08306B","Jan" = "#2171B5","Feb" = "#6BAED6","Mar" = "#74C476","Apr" = "#31A354","May" = "#006D2C",
                 "Jun" = "#FB6A4A","Jul" = "#DE2D26","Aug" = "#A50F15", "Sep" = "#DFC27D","Oct" = "#BF812D","Nov" = "#8C510A")
month.lty = rep(c(1,2,4), times = 4)
names(month.lty) = month.abb

ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.sparse, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 29),
        strip.text   = element_text(size = 25),
        legend.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25)) +
  
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("weather temperature/figures/sd_sparse_",data.example[[k]],".png"), width = 28, height = 18, units = "cm", dpi = 300)

ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.dense, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  guides(colour      = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  
  theme(plot.title   = element_text(size =25),
        legend.text  = element_text(size =29,  color = "transparent"),
        strip.text   = element_text(size = 25, color = "transparent"),
        legend.title = element_text(size = 30, color = "transparent"),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        legend.position = "right",
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("weather temperature/figures/sd_dense_",data.example[[k]],".png"), width = 28, height = 18, units = "cm", dpi = 300)








#############################################################################################
##### Plot of long run standard deviation of each month #####################################
#############################################################################################


ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.sparse.lr, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  
  theme(plot.title   = element_text(size =25),
        legend.text  = element_text(size =29),
        strip.text   = element_text(size = 25),
        legend.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25)) +
  
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("weather temperature/figures/sd_sparse_long_run_",data.example[[k]],".png"), width = 28, height = 18, units = "cm", dpi = 300)

ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.dense.lr, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  guides(colour      = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 29, color = "transparent"),
        strip.text   = element_text(size = 25, color = "transparent"),
        legend.title = element_text(size = 30, color = "transparent"),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        legend.position = "right",
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("weather temperature/figures/sd_dense_long_run_",data.example[[k]],".png"), width = 28, height = 18, units = "cm", dpi = 300)










#######################################################################
#### Plot of dense and sparse mean function for each month ############
#######################################################################


# Loading bandwidths for mean and difference estimation
Bandwidths = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_",data.example[[k]],".rds"))

# Mean and difference function estimation
est = est.results(data.s.34h,data.d.34h,Bandwidths,from = 1,to = 12)

if(k == 1){
  
ggplot() +
  geom_line(aes(x = TIME, y = ESTIMATE), data =est$dense, color = "turquoise", lty = 1, size = 1.5) +
  geom_line(mapping = aes(x = TIME, y = ESTIMATE), data = est$sparse, color = "darkblue", size = 1.2, lty = 2, show.legend = F) +
  labs(x = "Time", y = "Temperature in °C",title = bquote(.(data.example.pic[[k]]) *" (Germany): Estimation of " *mu^{"[d]"} * " and " *mu^{"[s]"})) +
    
  theme(plot.title = element_text(size =17),
        legend.text = element_text(size =10),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 11),
        axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 12),     
        axis.text.y  = element_text(size = 12)) +
    
  ylim(-4,27)+
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))+
  facet_wrap(~ MONTH, ncol = 6, nrow = 2)
ggsave(paste0("weather temperature/figures/means_",data.example[[k]],".png"), width = 30, height = 14, units = "cm", dpi = 300)

}else{
  
ggplot() +
  geom_line(aes(x = TIME, y = ESTIMATE), data =est$dense, color = "turquoise", lty = 1, size = 1.5) +
  geom_line(mapping = aes(x = TIME, y = ESTIMATE), data = est$sparse, color = "darkblue", size = 1.2, lty = 2, show.legend = F) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("Estimation of " *mu^{"[d]"} * " and " *mu^{"[s]"}*phantom(integral(delta)))) +
    
  theme(plot.title = element_text(size = 22),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),     
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 17),     
        axis.text.y  = element_text(size = 17)) +
    
  ylim(-4,27)+
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("0", "10", "20"))+
  facet_wrap(~ MONTH, ncol = 12, nrow = 1)
ggsave(paste0("weather temperature/figures/means_",data.example[[k]],".png"), width = 46, height = 10, units = "cm", dpi = 300)

}





########################################################################
#### Plot of difference function and integral with confidence bands ####
########################################################################

integral.conf = readRDS(paste0("weather temperature/confidence bands/Results/CB_integral_",data.example[[k]],".rds"))
delta.conf    = readRDS(paste0("weather temperature/confidence bands/Results/CB_delta_",data.example[[k]],".rds"))


if(k == 1){
  
ggplot() +
  labs(x = "Time", y = "Temperature in °C",title = bquote(.(data.example.pic[[k]]) *" (Germany): Estimation of  "* delta * "  and  " * integral(delta) * " d" *lambda)) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = integral.conf, fill = "darkred", col = NA, alpha = 0.3,size = 0.2)+
  geom_line(aes(x = TIME, y = ESTIMATE, color = ESTIMATE), data =est$delta_int, lty = 2, size = 1.8, alpha = 1, show.legend = T) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = delta.conf, fill = "grey6", col = NA, alpha = 0.2,size = 0.5)+
  geom_line(mapping = aes(x = TIME, y = ESTIMATE, color = ESTIMATE), data = est$delta, size = 1.8, show.legend = F) +
  scale_color_gradient2(mid = "darkblue", high = "red", midpoint = -1, oob = scales::squish, limits = c(-1, 3), name = "°C") +
    
    theme(plot.title = element_text(size   = 15),
          legend.text = element_text(size  = 12),
          strip.text = element_text(size   = 14),
          legend.title = element_text(size = 13),
          axis.title.x = element_text(size = 13),     
          axis.title.y = element_text(size = 14),
          axis.text.x  = element_text(size = 11),     
          axis.text.y  = element_text(size = 11)) +
    
  scale_y_continuous(breaks = c(-1.5,0.,1.5,3,4.5), limits = c(-1.7, 5.9)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))+
  facet_wrap(~ MONTH, ncol = 6, nrow = 2)
ggsave(paste0("weather temperature/figures/difference_",data.example[[k]],".png"), width = 30, height = 14, units = "cm", dpi = 300)

}else{

ggplot() +
  labs(x = "Time", y = "Temperature in °C",title = bquote("Estimation of  "* delta * "  and  " * integral(delta) * " d" *lambda)) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = integral.conf, fill = "darkred", col = NA, alpha = 0.3,size = 0.2)+
  geom_line(aes(x = TIME, y = ESTIMATE, color = ESTIMATE), data =est$delta_int, lty = 2, size = 1.8, alpha = 1, show.legend = F) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = delta.conf, fill = "grey6", col = NA, alpha = 0.2,size = 0.5)+
  geom_line(mapping = aes(x = TIME, y = ESTIMATE, color = ESTIMATE), data = est$delta, size = 1.8, show.legend = F) +
  scale_color_gradient2(mid = "darkblue", high = "red", midpoint = -1, oob = scales::squish, limits = c(-1, 3), name = "°C") +
    
    theme(plot.title = element_text(  size = 22),
          strip.text = element_text(  size = 20),
          axis.title.x = element_text(size = 20),     
          axis.title.y = element_text(size = 20),
          axis.text.x  = element_text(size = 17),     
          axis.text.y  = element_text(size = 17)) +
    
  scale_y_continuous(breaks = c(0,2.5,5), limits = c(-1.7, 5.9)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("0", "10", "20"))+
  facet_wrap(~ MONTH, ncol = 12, nrow = 1)
ggsave(paste0("weather temperature/figures/difference_",data.example[[k]],".png"), width = 46, height = 10, units = "cm", dpi = 300)

}









#####################################################################
# Calculating diagonal of long run covariance matrix for each month #
#####################################################################

lr.Gamma.d.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                         function(j){if(j >= 2){
                                                                           (cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                         }else{cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                         }))
return(diag(lr.cov))})

lr.Gamma.s.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                         function(j){if(j >= 2){
                                                                           (cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]+t(cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                         }else{cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                         }))
return(diag(lr.cov))})

lr.Gamma = lapply(1:12, function(m) { sample.dense  = data.d.24h.hourly %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d.24h.hourly)[2])
                                      sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
                                      sample.sparse = data.s.24h %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s.24h)[2])
                                      sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
                                      
                                      return(lr.Gamma.s.month[[m]] + dim(sample.sparse)[1] / dim(sample.dense)[1] * lr.Gamma.d.month[[m]])
                                     })


centered.delta.conf = readRDS(paste0("weather temperature/confidence bands/Results/CB_centered_delta_",data.example[[k]],".rds"))

# Comparison: Confidence bands of difference function (not centered)
q.list              = readRDS(paste0("weather temperature/confidence bands/Results/quantile_",data.example[[k]],".rds"))
delta.conf.compare  = CB(data.s.34h, est, lr.Gamma, q.list[1,], center = T)


if(k == 1){

ggplot() +
  labs(x = "Time", y = "Temperature in °C",title = bquote(.(data.example.pic[[k]]) *" (Germany): Estimation of  "* delta - integral(delta) * " d" *lambda)) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = centered.delta.conf,  fill  = "grey6",  col = NA, alpha = 0.3,size = 0.5)+
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = delta.conf.compare, fill  = "grey6",  col = NA, alpha = 0.2,size = 0.5, lty = 2)+
  geom_line(mapping = aes(x = TIME, y = ESTIMATE), data = centered.delta.conf, color = "red", size = 0.8, show.legend = F, linetype = 2) +
  
    theme(plot.title = element_text(size =17),
          legend.text = element_text(size =10),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 11),
          axis.title.x = element_text(size = 14),     
          axis.title.y = element_text(size = 15),
          axis.text.x  = element_text(size = 12),     
          axis.text.y  = element_text(size = 12)) +
  
  ylim(-3.6, 4)+
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))+
  facet_wrap(~ MONTH, ncol = 6, nrow = 2)
ggsave(paste0("weather temperature/figures/centered_difference_",data.example[[k]],".png"), width = 30, height = 14, units = "cm", dpi = 300)

}else{

ggplot() +
  labs(x = "Time", y = "Temperature in °C",title = bquote("Estimation of  "* delta - integral(delta) * " d" *lambda)) +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = centered.delta.conf, fill  = "grey6",  col = NA, alpha = 0.3,size = 0.5)+
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = delta.conf.compare,  fill  = "grey6",  col = NA, alpha = 0.2,size = 0.5, lty = 2)+
  geom_line(mapping = aes(x = TIME, y = ESTIMATE), data = centered.delta.conf, color = "red", size = 0.9, show.legend = F, linetype = 2) +
  
    theme(plot.title = element_text(size = 22),
          strip.text = element_text(size = 20),
          axis.title.x = element_text(size = 20),     
          axis.title.y = element_text(size = 20),
          axis.text.x  = element_text(size = 17),     
          axis.text.y  = element_text(size = 17)) +
  
  scale_y_continuous(breaks = c(-3,0,3), limits = c(-3.6, 4)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("0", "10", "20"))+
  facet_wrap(~ MONTH, ncol = 12, nrow = 1)
ggsave(paste0("weather temperature/figures/centered_difference_",data.example[[k]],".png"), width = 46, height = 10, units = "cm", dpi = 300)

}






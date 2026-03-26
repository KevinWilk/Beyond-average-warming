library(ggplot2)
library(reshape2)
library(locpol)
library(biLocPol)
library(interp)
library(stats)
library(future)
library(future.apply)
library(parallel)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tibble)
library(plotly)
library(gridExtra)
library(CompQuadForm)
library(tensorA)



l = 1

data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")

load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))
#load(paste0("Application/Data examples/function_to_load.RData"))


########################################################
# Estimating lag covariance kernels for lag = 0,...,12 #
########################################################


cov.s.list = list()
cov.d.list = list()


options(future.globals.maxSize = 20 * 1024^3) 
plan(multisession, workers = future::availableCores()-2)


for(m in 1:12){
  
    sample.dense  = data.d.34h |> 
                      filter(MONTH == month.name[m]) |> 
                        dplyr::select(1,4:dim(data.d.34h)[2])
    
    sample.dense = sample.dense[rowSums(is.na(sample.dense)) == 0,] 
    pd = dim(sample.dense)[2] - 1
  
    sample.sparse = data.s.34h |> 
                      filter(MONTH == month.name[m]) |> 
                        dplyr::select(1,4:dim(data.s.34h)[2])
    
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]  
    p = dim(sample.sparse)[2] - 1
  
    lag.k.Gamma   = list()
    lag.k.Gamma.d = list()  
  
    file.s = paste0("weather temperature/kernel weights/full_w_s_lag0.rds")
    file.d = paste0("weather temperature/kernel weights/full_w_d_lag0," ",.rds")  ##HIER PROBLEME
  
    w   = readRDS(file.s)
    wd  = readRDS(file.d)
  
    lag.Gamma = eval.weights(w,observation.transformation(sample.sparse[,-1],grid.type = "less", periodic = T, m = 24))
    lag.k.Gamma[[1]] = lag.Gamma

    lag.Gamma.d   = eval.weights(wd,observation.transformation(sample.dense[,-1],grid.type = "less", periodic = T, m = 144))
    lag.k.Gamma.d[[1]] = lag.Gamma.d
  
    rm(w,lag.Gamma)
    rm(wd,lag.Gamma.d)
  
    lag.k.Gamma.part2 = future_lapply(1:2, function(k) {
      file.s = paste0("weather temperature/kernel weights/full_w_s_lag",k,".rds")
      w.lag = readRDS(file.s)
      n.year = unique(sample.sparse$Year) 
      lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(sample.sparse[which(sample.sparse$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 24), lag = k)}))/length(n.year)
      lag.Gamma
    })
    lag.k.Gamma = c(lag.k.Gamma,lag.k.Gamma.part2)
  
    #lag.k.Gamma.d.part2 = future_lapply(1:2, function(k) {  
    #  file.d = paste0("weather temperature R files/kernel weights/",data.example[[l]],month.name[m],"/dense/full_w.lag",k,".rds")
    #  w.lag = readRDS(file.d)
    #  n.year = unique(sample.dense$Year) 
    #  lag.Gamma = Reduce(`+`,lapply(n.year,  function(j){eval.weights(w.lag, observation.transformation(sample.dense[which(sample.dense$Year %in% j), -1], lag = k, grid.type = "lesseq", periodic = T, m = 144), lag = k)}))/length(n.year)
    #  lag.Gamma
    #})
    #lag.k.Gamma.d = c(lag.k.Gamma.d,lag.k.Gamma.d.part2)
  
    cov.s.list[[m]]   = lag.k.Gamma
    #cov.d.list[[m]]   = lag.k.Gamma.d
  
    rm(lag.k.Gamma,p)
    rm(lag.k.Gamma.d,pd)
  
    rm(sample.dense,sample.sparse,lag.k.Gamma.d.part2,lag.k.Gamma.part2)
    print(paste0("done: ", month.name[m]))
}

plan(sequential)







#######################
# Finding maximal lag #
#######################

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


 
options(future.globals.maxSize = 20 * 1024^3) #laptop hat 21,9 GB zur Verfügung
plan(multisession, workers = future::availableCores()-2)
test.s.val = lapply(1:12, function(month) {inference.autocovariance.test(data.s.24h   |> select(1:27) , m = month, alpha = 0.95, max.lag = 12)})
test.d.val = lapply(1:12, function(month) {inference.autocovariance.test(dense.hourly |> select(1:27) , m = month, alpha = 0.95, max.lag = 12)}) 
plan(sequential)
test.s.val = data.frame(lag = rep(c(0.6,2:11,12.5), times = 12), value = unlist(test.s.val), MONTH = factor(rep(month.name[1:length(test.s.val)], lengths(test.s.val)), levels = month.name))
test.d.val = data.frame(lag = rep(c(0.6,2:11,12.5), times = 12), value = unlist(test.d.val), MONTH = factor(rep(month.name[1:length(test.d.val)], lengths(test.d.val)), levels = month.name))

ggplot() +
  geom_bar(aes(x = lag,  y     = value), stat = "identity", fill = "turquoise",width = 0.4, data = rho.hat.d.df, alpha = 0.7) +
  geom_line(aes(x = lag, y = value), test.d.val, col = "green3", size = 1.1, lty = 2)+
  geom_line(aes(x = lag+0.3, y = value), test.s.val, col = "blue", size = 1.1, lty = 2)+
  geom_bar(aes(x = lag+0.3, y = value),  stat = "identity", fill = "darkblue", width = 0.4, data = rho.hat.s.df, alpha = 0.5) +
  labs(x = "Lag", y = bquote(hat(rho)[Lag]) ,title = "") +
  theme(plot.title = element_text(size =15),
        legend.text = element_text(size =15),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 15),     
        axis.text.y  = element_text(size = 15))+
  scale_y_continuous(breaks = c(0,0.25,0.5), limits = c(0,0.7))+
  scale_x_continuous(breaks = c(1,4,8,12), limits = c(0.5,12.8))+
  facet_wrap(~ MONTH)
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/lag test.png"), width = 32, height = 18, units = "cm", dpi = 300)

month.max.lag

month.d.max.lag = c(2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 2, 1);month.s.max.lag = c(2, 1, 2, 2, 3, 4, 2, 2, 2, 2, 1, 1)#Berlin
month.d.max.lag = c(2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2);month.s.max.lag = c(2, 0, 2, 4, 2, 4, 4, 2, 0, 1, 1, 1)#Frankfurt
month.d.max.lag = c(2, 2, 1, 1, 1, 2, 1, 2, 2, 1, 1, 1);month.s.max.lag = c(2, 1, 1, 2, 1, 4, 1, 2, 0, 1, 1, 1)#Hamburg
month.d.max.lag = c(2, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1);month.s.max.lag = c(2, 1, 2, 2, 3, 5, 2, 2, 1, 1, 1, 1)#Leipzig
month.d.max.lag = c(2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2);month.s.max.lag = c(3, 1, 2, 3, 3, 5, 2, 2, 1, 2, 1, 2)#Muenchen
month.max.lag = apply(rbind(month.d.max.lag,month.s.max.lag),2,max) 


#############################################################################################
##### Plot of sd of each month ##############################################################
#############################################################################################

start.24 = 21; end.24 = 117 

sqrt.d.month = lapply(1:12, function(m) {sqrt(diag(cov.d.list[[m]][[1]][start.24:end.24,start.24:end.24]))})
sqrt.s.month = lapply(1:12, function(m) {sqrt(diag(cov.s.list[[m]][[1]][start.24:end.24,start.24:end.24]))})
lr.sqrt.d.month = lapply(1:12, function(m) {sqrt(Reduce(`+`,lapply(1:(month.max.lag[m]+1),  function(j){diag(cov.d.list[[m]][[j]][start.24:end.24,start.24:end.24]) * ifelse(j >= 2, 2*(1-(j-1)/(month.max.lag[m]+1)), 1) })))})
lr.sqrt.s.month = lapply(1:12, function(m) {sqrt(Reduce(`+`,lapply(1:(month.max.lag[m]+1),  function(j){diag(cov.s.list[[m]][[j]][start.24:end.24,start.24:end.24]) * ifelse(j >= 2, 2*(1-(j-1)/(month.max.lag[m]+1)), 1) })))})
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
  theme(plot.title = element_text(size =25),
        legend.text = element_text(size =29),
        strip.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25))+
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd sparse monthly.png"), width = 28, height = 18, units = "cm", dpi = 300)

ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.dense, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  guides(colour   = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  theme(plot.title = element_text(size =25),
        legend.text = element_text(size =29,color = "transparent"),
        strip.text = element_text(size = 25,  color = "transparent"),
        legend.title = element_text(size = 30,color = "transparent"),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        legend.position = "right",
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = "transparent", color = NA))+
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd dense monthly.png"), width = 28, height = 18, units = "cm", dpi = 300)


ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.sparse.lr, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  theme(plot.title = element_text(size =25),
        legend.text = element_text(size =29),
        strip.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25))+
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd sparse monthly long run.png"), width = 28, height = 18, units = "cm", dpi = 300)

ggplot() +
  geom_line(mapping = aes(x = TIME, y = SD.dense.lr, col = MONTH, linetype = MONTH), data = sqrt.month, size = 2, show.legend = T, alpha = 0.9) +
  labs(x = "Time", y = "Temperature in °C",title = bquote("")) +
  scale_color_manual(values = month.colors, name = "Month") +
  scale_linetype_manual(values = month.lty,limits = month.abb, name = "Month") +
  guides(colour   = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  theme(plot.title = element_text(size =25),
        legend.text = element_text(size =29,color = "transparent"),
        strip.text = element_text(size = 25,  color = "transparent"),
        legend.title = element_text(size = 30,color = "transparent"),
        axis.title.x = element_text(size = 25),     
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 25),     
        axis.text.y  = element_text(size = 25),
        legend.position = "right",
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = "transparent", color = NA))+
  scale_y_continuous(breaks = c(0, 2.5, 5 , 7.5, 10), limits = c(0, 9.5)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd dense monthly long run.png"), width = 28, height = 18, units = "cm", dpi = 300)



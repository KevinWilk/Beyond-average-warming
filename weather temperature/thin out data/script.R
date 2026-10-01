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


bw.s_month = function(data.sparse,data.dense, h.seq){
  
  bw_list = matrix(NA, nrow = 3, ncol = 12,dimnames = list(c("dense", "delta", "sparse"), month.name)) |> as.data.frame()
  
  results = future_sapply(c(1,4,7,10), function(m){
    
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    
    sample.dense  = data.dense %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data.dense)[2])
    sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
    K.dense  = length(unique(sample.dense$Year))
    
    sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m.neighbor])|> dplyr::select(1,4:dim(data.sparse)[2])
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
    K.sparse = length(unique(sample.sparse$Year))
    
    bw.dense      = round(readRDS(paste0("weather temperature/bandwidth selection/Results/bw_",data.example[[k]],".rds"))[1,m] * (34/24) , digits = 2)
    
    bw.delta      = k.fold.hc.cv(sample = sample.sparse, h.seq = h.seq, deg = 2, diff = T, K = K.sparse, sample.dense, bw.dense) 
    bw.sparse     = k.fold.hc.cv(sample = sample.sparse, h.seq = h.seq, deg = 2, K = K.sparse)  
    

    c(bw.dense,bw.delta,bw.sparse)
    
  }, future.seed = T)
  
  bw_list[1, c(12,1:11)] = rep(results[1,],each = 3)
  bw_list[2, c(12,1:11)] = rep(results[2,],each = 3)
  bw_list[3, c(12,1:11)] = rep(results[3,],each = 3)
  
  return(bw_list)
}






options(future.globals.maxSize = 20 * 1024^3)                
plan(multisession, workers = 20) 

for(k in 1:4){
                                             
  data.example = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")      
  load(paste0("weather temperature/data sets/",data.example[[k]],".RData")) 

  Bandwidths = bw.s_month(data.s.24h[c(1:3,seq(4,28,2))], data.d.24h, h.seq = seq(0.18, 0.24, 0.01))
  file = paste0("weather temperature/thin out data/bw.2h_",data.example[[k]],".rds")
  saveRDS(Bandwidths, file)
  
  Bandwidths = bw.s_month(data.s.24h[c(1:3,seq(4,28,3))], data.d.24h, h.seq = seq(0.27, 0.34, 0.01))
  file = paste0("weather temperature/thin out data/bw.3h_",data.example[[k]],".rds")
  saveRDS(Bandwidths, file)
  
  print("save done")
  
  
}

plan(sequential)







###############################################################################
k = 1  # Change to: 1 (Berlin)                                               ##
#                   2 (Frankfurt am Main)                                    ##
#                   3 (Hamburg)                                              ##
#                   4 (Munich)                                               ##
data.example     = list("Berlin", "Frankfurt_Main", "Hamburg", "Munich")     ##
data.example.pic = list("Berlin", "Frankfurt am Main", "Hamburg", "Munich")  ##
load(paste0("weather temperature/data sets/",data.example[[k]],".RData"))    ##
###############################################################################



Bandwidths = readRDS(paste0("weather temperature/bandwidth selection/Results/bw_",data.example[[k]],".rds"))
est = est.results(data.s.34h, data.d.34h, Bandwidths, from = 1,to = 12)

Bandwidths.2h = readRDS(paste0("weather temperature/thin out data/bw.2h_",data.example[[k]],".rds"))
est.2h        = est.results(data.s.24h[c(1:3,seq(4,28,2))], data.d.24h, Bandwidths.2h, from = 1,to = 12, eval.grid = 0:96/96)

Bandwidths.3h = readRDS(paste0("weather temperature/thin out data/bw.3h_",data.example[[k]],".rds"))
est.3h        = est.results(data.s.24h[c(1:3,seq(4,28,3))], data.d.24h, Bandwidths.3h, from = 1,to = 12, eval.grid = 0:96/96)

centered.delta.conf = readRDS(paste0("weather temperature/confidence bands/Results/CB_centered_delta_",data.example[[k]],".rds"))



labs.legend = c(
  "res" = expression(hat(Delta)["145,25"]),
  "25"  = expression(bar(Delta)["145,25"]),
  "13"  = expression(bar(Delta)["145,13"]),
  "9"   = expression(bar(Delta)["145,9"])
)


ggplot() +
  labs(x = "Time", y = "Temperature in °C") +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = centered.delta.conf,  fill  = "grey6",  col = NA, alpha = 0.3,size = 0.5)+
  
  geom_line(aes(x = TIME, y = centered, colour = "25" , linetype = "25" ), data = est$delta_compare,    size = 1, show.legend = T) +
  geom_line(aes(x = TIME, y = centered, colour = "13" , linetype = "13" ), data = est.2h$delta_compare, size = 1, show.legend = T) +
  geom_line(aes(x = TIME, y = centered, colour = "9"  , linetype = "9"  ), data = est.3h$delta_compare, size = 1, show.legend = T) +
  geom_line(aes(x = TIME, y = centered, colour = "res", linetype = "res"), data = est$delta,            size = 1, show.legend = T) +
  
  scale_colour_manual(expression("Estimators:"), 
                    values = c("res" = "red","25"  = "blue", "13" = "green", "9" = "purple"),
                    labels = labs.legend,
                    breaks = c("res","25", "13", "9")) +
  
  scale_linetype_manual(expression("Estimators:"),
                        values = c( "res"  = 2, "25" = 1, "13" = 3, "9"  = 4),
                        labels = labs.legend,
                        breaks = c("res","25", "13", "9")) +
  
  
  theme(plot.title = element_text(size =17),
        legend.text = element_text(size =16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 12),     
        axis.text.y  = element_text(size = 12),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width  = unit(1.1, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  
  scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.5, 1.6)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))+
  facet_wrap(~ MONTH, ncol = 6, nrow = 2)

ggsave(paste0("weather temperature/thin out data/thin_out_",data.example[[k]],".png"), width = 35, height = 11, units = "cm", dpi = 300)



count.h  = 0
count.2h = 0
count.3h = 0

for(m in 1:12){
  
  compare = est$delta_compare   |> filter(MONTH == month.name[m]) |> select(centered)
  UP      = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(UP)
  LOW     = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(LO)
  p25 = bind_cols(compare, LOW, UP) |> mutate(between = LO <= centered & centered <= UP)
  if(any(p25$between == F)){count.h = count.h + 1}
  
  compare = est.2h$delta_compare   |> filter(MONTH == month.name[m]) |> select(centered)
  UP      = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(UP)
  LOW     = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(LO)
  p13 = bind_cols(compare, LOW, UP) |> mutate(between = LO <= centered & centered <= UP)
  if(any(p13$between == F)){count.2h = count.2h + 1}
  
  compare = est.3h$delta_compare   |> filter(MONTH == month.name[m]) |> select(centered)
  UP      = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(UP)
  LOW     = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(LO)
  p9 = bind_cols(compare, LOW, UP) |> mutate(between = LO <= centered & centered <= UP)
  if(any(p9$between == F)){count.3h = count.3h + 1}
  
  
}

count.h; count.2h; count.3h


ggplot() +
  labs(x = "Time", y = "Temperature in °C") +
  geom_ribbon(aes(x = TIME, ymin = LO, ymax = UP), data = centered.delta.conf,  fill  = "grey6",  col = NA, alpha = 0.3,size = 0.5)+
  
  geom_line(aes(x = TIME, y = centered, colour = "25", linetype = "25"), data = est$delta,    size = 1, show.legend = T) +
  geom_line(aes(x = TIME, y = centered, colour = "13", linetype = "13"), data = est.2h$delta, size = 1, show.legend = T) +
  geom_line(aes(x = TIME, y = centered, colour = "9",  linetype = "9" ), data = est.3h$delta, size = 1, show.legend = T) +

  scale_colour_manual(expression("Estimators:"), 
                      values = c("25" = "red", "13" = "green", "9" = "purple"),
                      labels = c("25"  = expression(hat(Delta)["145,25"]),
                                 "13"  = expression(hat(Delta)["145,13"]),
                                 "9"   = expression(hat(Delta)["145,9"])),
                      breaks = c("25","13","9")) +
  
  scale_linetype_manual(expression("Estimators:"),
                        values = c("25" = 2, "13" = 3, "9"  = 4),
                        labels = c("25" = expression(hat(Delta)["145,25"]),
                        "13" = expression(hat(Delta)["145,13"]),
                        "9"  = expression(hat(Delta)["145,9"])),
                        breaks = c("25", "13", "9")) +
  
  theme(plot.title = element_text(size =17),
        legend.text = element_text(size =16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),     
        axis.title.y = element_text(size = 15),
        axis.text.x  = element_text(size = 12),     
        axis.text.y  = element_text(size = 12),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width  = unit(1.1, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  
  scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.5, 1.6)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00"))+
  facet_wrap(~ MONTH, ncol = 6, nrow = 2)

ggsave(paste0("weather temperature/thin out data/thin_out_residuals_",data.example[[k]],".png"), width = 35, height = 11, units = "cm", dpi = 300)


count.2h = 0
count.3h = 0

for(m in 1:12){
  
  compare = est.2h$delta        |> filter(MONTH == month.name[m]) |> select(centered)
  UP      = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(UP)
  LOW     = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(LO)
  p13 = bind_cols(compare, LOW, UP) |> mutate(between = LO <= centered & centered <= UP)
  if(any(p13$between == F)){count.2h = count.2h + 1}
  
  compare = est.3h$delta        |> filter(MONTH == month.name[m]) |> select(centered)
  UP      = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(UP)
  LOW     = centered.delta.conf |> filter(MONTH == month.name[m]) |> select(LO)
  p9 = bind_cols(compare, LOW, UP) |> mutate(between = LO <= centered & centered <= UP)
  if(any(p9$between == F)){count.3h = count.3h + 1}
  
  
}

count.2h; count.3h


####################################################################################################
####################################################################################################
###############              #######################################################################
###############    Figures   #######################################################################
###############              #######################################################################
####################################################################################################
####################################################################################################



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



#############################################################################################
##### Plot of inference autocovariance test results of each month ###########################
#############################################################################################

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

ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/lag test.png"), width = 32, height = 18, units = "cm", dpi = 300)



#############################################################################################
##### Plot of standard deviation (lag = 0) of each month ####################################
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
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd sparse monthly.png"), width = 28, height = 18, units = "cm", dpi = 300)

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
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd dense monthly.png"), width = 28, height = 18, units = "cm", dpi = 300)



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
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd sparse monthly long run.png"), width = 28, height = 18, units = "cm", dpi = 300)

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
ggsave(paste0("Application/Data examples/",data.example[[l]],"/pictures/sd dense monthly long run.png"), width = 28, height = 18, units = "cm", dpi = 300)

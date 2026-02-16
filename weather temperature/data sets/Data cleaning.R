
# ---- Setup / Packages ----
# renv (R-enviroment) is managed outside this script:
#   - run renv::restore() after cloning


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
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)                    
library(rdwd)


 
for(k in 1:4){
  
  cities = data.frame(
    city = c("Berlin","Frankfurt_Main","Hamburg","Munich"), 
    sparse.station = c("00433;00424", "01425;01421" , "01975" , "03390"), 
    dense.station  = c("00433", "01420" , "01975" , "03379") 
  )
  
  ####################################
  ####################################
  ####                            ####
  ####  Generating dense data set ####
  ####                            ####
  ####################################
  ####################################
  
  link = selectDWD(id = cities$dense.station[k], res = "10_minutes", var = "air_temperature", per = "historical")

  data.d.list  = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
  data.d.part1 = data.d.list[[1]]
  data.d.part2 = subset(data.d.list[[2]], select = -eor)
  data.d.part3 = subset(data.d.list[[3]], select = -eor)
  data.d.part4 = subset(data.d.list[[4]], select = -eor)
  data.d.part4 = data.d.part4[data.d.part4$MESS_DATUM <= as.POSIXct("2024-12-31 23:50:00", tz = "UTC"), ]

  link = selectDWD(id = cities$dense.station[k], res = "10_minutes", var = "air_temperature", per = "recent")

  data.d.part5 = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
  data.d.part5 = subset(data.d.part5[data.d.part5$MESS_DATUM >= as.POSIXct("2025-01-01 00:00:00", tz = "UTC") & 
                                       data.d.part5$MESS_DATUM <= as.POSIXct("2026-01-01 23:50:00", tz = "UTC"), ], select = -eor)

  data.d =  rbind(data.d.part1,data.d.part2,data.d.part3,data.d.part4,data.d.part5)
  rm(link, data.d.list, data.d.part1, data.d.part2, data.d.part3, data.d.part4, data.d.part5)

  data.d$Year  = year(data.d$MESS_DATUM)
  data.d$MONTH = month(data.d$MESS_DATUM)
  data.d$MONTH = factor(data.d$MONTH, levels = 1:12, labels = month.name)
  data.d$DAY   = day(data.d$MESS_DATUM)
  data.d$TIME  = as_hms(data.d$MESS_DATUM)

  data.d = data.d[data.d$MESS_DATUM >= as.POSIXct("1999-12-31 00:00:00", tz = "UTC"), ]
  
  data.d = data.d |> 
              mutate(TIME = as.character(TIME)) |> 
                dplyr::select(Year, MONTH, DAY, TIME, TT_10) |> 
                  pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_10,values_fn = first)

  data.d.before.5h           = data.d[-c(dim(data.d)[1]-1,dim(data.d)[1]),c(118:147)] 
  colnames(data.d.before.5h) = paste("Before",colnames(data.d.before.5h))
  data.d.after.5h            = data.d[-c(1,2),c(4:34)]    
  colnames(data.d.after.5h)  = paste("After",colnames(data.d.after.5h))
  
  
  
  ######################################################################
  # dense data set from 19:00 to 05:00 (205 design points, length 34h) #
  ######################################################################
  
  data.d.34h                 = cbind(data.d[-c(1,dim(data.d)[1]),1:3],
                                     data.d.before.5h,
                                        data.d[-c(1,dim(data.d)[1]),-c(1:3)],
                                            data.d.after.5h)
  
  
  
  rm(data.d.before.5h, data.d.after.5h)

  data.d.before.3h           = data.d[-c(dim(data.d)[1]-1,dim(data.d)[1]),c(130:147)] 
  colnames(data.d.before.3h) = paste("Before",colnames(data.d.before.3h))
  data.d.after.3h            = data.d[-c(1,2),c(4:22)]     
  colnames(data.d.after.3h)  = paste("After",colnames(data.d.after.3h))
  
  
  
  ##########################################################################################
  # dense data set from 21:00 to 03:00 (181 design points, length 30h)                     #
  # Only usage: Bandwidth selection of the lag 0 covariance kernels                        #
  # Reason: 34-hour dataset (205 design points) exceeds available computational resources  #
  ##########################################################################################
  
  data.d.30h                 = cbind(data.d[-c(1,dim(data.d)[1]),1:3],
                                     data.d.before.3h,
                                        data.d[-c(1,dim(data.d)[1]),-c(1:3)],
                                            data.d.after.3h)
  
  
  
  
  rm(data.d.before.3h, data.d.after.3h)
  
  keep.cols  = names(data.d.34h)[grepl("^(Year|MONTH|DAY)$", names(data.d.34h)) | 
                                   names(data.d.34h) == "After 00:00:00" |
                                      (!grepl("^Before ", names(data.d.34h)) & !grepl("^After ", names(data.d.34h)))]
  
  

  
  
  ####################################################################################################################
  # dense data set hourly from 00:00 to 00:00 (25 design points, length 24h)                                         #
  # Only usage: Test to assess the cumulative significance of empirical lagged covariance kernel                     #
  # Note: Tensor products of size 145 are computationally too demanding; therefore, hourly aggregated data were used #          
  ####################################################################################################################
  
  data.d.24h = data.d.30h |> select(all_of(keep.cols))
  
  hour.agg = function(data, pattern = NULL) {
                data |>
                    pivot_longer(cols = -c(Year, MONTH, DAY),names_to = "name",values_to = "temp") %>%
                      { if (!is.null(pattern)) filter(., str_detect(name, pattern)) else . } %>%
                        mutate(time = str_extract(name, "\\d{2}:\\d{2}:\\d{2}$")) %>%
                          filter(!is.na(time)) |>mutate(hour = as.integer(str_sub(time, 1, 2)), hourly = sprintf("%02d:00:00", hour)) %>%
                            group_by(Year, MONTH, DAY, hourly) %>%
                              summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
                                pivot_wider(names_from = hourly, values_from = temp) %>%
                                  arrange(Year, MONTH, DAY)
  }
  
  data.d.24h.hourly = hour.agg(data.d.24h, pattern = "^(?!Before\\s|After\\s).*\\d{2}:\\d{2}:\\d{2}$") |>
                        mutate(`After 00:00:00` = `00:00:00`)  
  
  rm(data.d.24h, hour.agg)
  
  
  
  
  #####################################
  #####################################
  ####                             ####
  ####  Generating sparse data set ####
  ####                             ####
  #####################################
  #####################################
  
  if(cities$city[k] == "Berlin"){
    
    link = selectDWD(id = strsplit(cities$sparse.station[k], ";")[[1]][1] , res = "hourly", var = "air_temperature", per = "historical") 

    data.s.part1 = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
    rm(link)
    
    data.s.part1$Year  = year(data.s.part1$MESS_DATUM)
    data.s.part1$MONTH = month(data.s.part1$MESS_DATUM)
    data.s.part1$MONTH = factor(data.s.part1$MONTH, levels = 1:12, labels = month.name)
    data.s.part1$DAY   = day(data.s.part1$MESS_DATUM)
    data.s.part1$TIME  = as_hms(data.s.part1$MESS_DATUM)
    
    data.s.part1 = data.s.part1[data.s.part1$MESS_DATUM >= as.POSIXct("1951-12-31 00:00:00", tz = "UTC") & 
                                  data.s.part1$MESS_DATUM <= as.POSIXct("1970-12-31 23:00:00", tz = "UTC"), ]
    
    data.s.part1 = data.s.part1 |> 
                    mutate(TIME = as.character(TIME)) |> 
                      dplyr::select(Year, MONTH, DAY, TIME, TT_TU) |> 
                        pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_TU,values_fn = first)
    
    
    link = selectDWD(id = strsplit(cities$sparse.station[k], ";")[[1]][2] , res = "hourly", var = "air_temperature", per = "historical") 
    
    data.s.part2 = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
    rm(link)
    
    data.s.part2$Year  = year(data.s.part2$MESS_DATUM)
    data.s.part2$MONTH = month(data.s.part2$MESS_DATUM)
    data.s.part2$MONTH = factor(data.s.part2$MONTH, levels = 1:12, labels = month.name)
    data.s.part2$DAY   = day(data.s.part2$MESS_DATUM)
    data.s.part2$TIME  = as_hms(data.s.part2$MESS_DATUM)
    
    data.s.part2 = data.s.part2[data.s.part2$MESS_DATUM >= as.POSIXct("1971-01-01 00:00:00", tz = "UTC") & 
                                  data.s.part2$MESS_DATUM <= as.POSIXct("1973-01-01 23:00:00", tz = "UTC"), ]
    
    data.s.part2 = data.s.part2 |> 
      mutate(TIME = as.character(TIME)) |> 
      dplyr::select(Year, MONTH, DAY, TIME, TT_TU) |> 
      pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_TU,values_fn = first)
    
    data.s = rbind(data.s.part1,data.s.part2)
    rm(data.s.part1,data.s.part2)
    
  }else if(cities$city[k] == "Frankfurt_Main"){
  
    link = selectDWD(id = strsplit(cities$sparse.station[k], ";")[[1]][1] , res = "hourly", var = "air_temperature", per = "historical") 
    
    data.s.part1 = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
    rm(link)
    
    data.s.part1$Year  = year(data.s.part1$MESS_DATUM)
    data.s.part1$MONTH = month(data.s.part1$MESS_DATUM)
    data.s.part1$MONTH = factor(data.s.part1$MONTH, levels = 1:12, labels = month.name)
    data.s.part1$DAY   = day(data.s.part1$MESS_DATUM)
    data.s.part1$TIME  = as_hms(data.s.part1$MESS_DATUM)
    
    data.s.part1 = data.s.part1[data.s.part1$MESS_DATUM >= as.POSIXct("1951-12-31 00:00:00", tz = "UTC") & 
                                  data.s.part1$MESS_DATUM <= as.POSIXct("1961-12-31 23:00:00", tz = "UTC"), ]
    
    data.s.part1 = data.s.part1 |> 
      mutate(TIME = as.character(TIME)) |> 
      dplyr::select(Year, MONTH, DAY, TIME, TT_TU) |> 
      pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_TU,values_fn = first)
    
    
    link = selectDWD(id = strsplit(cities$sparse.station[k], ";")[[1]][2] , res = "hourly", var = "air_temperature", per = "historical") 
    
    data.s.part2 = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
    rm(link)
    
    data.s.part2$Year  = year(data.s.part2$MESS_DATUM)
    data.s.part2$MONTH = month(data.s.part2$MESS_DATUM)
    data.s.part2$MONTH = factor(data.s.part2$MONTH, levels = 1:12, labels = month.name)
    data.s.part2$DAY   = day(data.s.part2$MESS_DATUM)
    data.s.part2$TIME  = as_hms(data.s.part2$MESS_DATUM)
    
    data.s.part2 = data.s.part2[data.s.part2$MESS_DATUM >= as.POSIXct("1962-01-01 00:00:00", tz = "UTC") & 
                                  data.s.part2$MESS_DATUM <= as.POSIXct("1973-01-01 23:00:00", tz = "UTC"), ]
    
    data.s.part2 = data.s.part2 |> 
      mutate(TIME = as.character(TIME)) |> 
      dplyr::select(Year, MONTH, DAY, TIME, TT_TU) |> 
      pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_TU,values_fn = first)
    
    
    data.s = rbind(data.s.part1,data.s.part2)
    rm(data.s.part1,data.s.part2)
    
  }else{
    
    link = selectDWD(id = cities$sparse.station[k], res = "hourly", var = "air_temperature", per = "historical") 
    
    data.s = dataDWD(link, force = FALSE, read = TRUE, quiet = TRUE)
    rm(link)
    
    data.s$Year  = year(data.s$MESS_DATUM)
    data.s$MONTH = month(data.s$MESS_DATUM)
    data.s$MONTH = factor(data.s$MONTH, levels = 1:12, labels = month.name)
    data.s$DAY   = day(data.s$MESS_DATUM)
    data.s$TIME  = as_hms(data.s$MESS_DATUM)
    
    data.s = data.s[data.s$MESS_DATUM >= as.POSIXct("1951-12-31 00:00:00", tz = "UTC") & 
                      data.s$MESS_DATUM <= as.POSIXct("1973-01-01 23:00:00", tz = "UTC"), ]
    
    data.s = data.s |> 
      mutate(TIME = as.character(TIME)) |> 
      dplyr::select(Year, MONTH, DAY, TIME, TT_TU) |> 
      pivot_wider(id_cols = c(Year, MONTH, DAY),names_from = TIME,values_from = TT_TU,values_fn = first)
    
  }
  
  
  data.s.before.5h           = data.s[-c(dim(data.s)[1]-1,dim(data.s)[1]),23:27]
  colnames(data.s.before.5h) = paste("Before",colnames(data.s.before.5h))
  data.s.after.5h            = data.s[-c(1,2),c(4:9)]                                               
  colnames(data.s.after.5h)  = paste("After",colnames(data.s.after.5h))
  
  
  
  ######################################################################
  # sparse data set from 19:00 to 05:00 (35 design points, length 34h) #
  ######################################################################
  
  data.s.34h                 = cbind(data.s[-c(1,dim(data.s)[1]),1:3],
                                      data.s.before.5h,
                                          data.s[-c(1,dim(data.s)[1]),-c(1:3)],
                                            data.s.after.5h)
  
  
  
  rm(data.s.before.5h, data.s.after.5h)
  
  
  
  
  ################################################################################################
  # sparse data set from 00:00 to 00:00 (25 design points, length 24h)                           #
  # Only usage: Test to assess the cumulative significance of empirical lagged covariance kernel #
  ################################################################################################
  
  keep.cols  = names(data.s.34h)[grepl("^(Year|MONTH|DAY)$", names(data.s.34h)) | 
                                   names(data.s.34h) == "After 00:00:00" |
                                   (!grepl("^Before ", names(data.s.34h)) & !grepl("^After ", names(data.s.34h)))]
  
  data.s.24h = data.s.34h |> select(all_of(keep.cols))
  
  
  
  
  
  filename = paste0("weather temperature/data sets/", cities$city[k], ".RData")
  rm(keep.cols, data.d, data.s, k, cities)
  
  save.image(file = filename)
  
}


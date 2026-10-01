
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
library(imputeTS)

 
for(k in 1:4){
  
  cities = data.frame(
    city = c("Berlin","Frankfurt_Main","Hamburg","Munich"), 
    sparse.station = c("00433;00424", "01425;01421" , "01975" , "03390"), 
    dense.station  = c("00433", "01420" , "01975" , "03379") 
  )
  #stations_id von_datum bis_datum Stationshoehe geoBreite geoLaenge   Stationsname                   Bundesland  Abgabe
  #00424:      19610101  19810101    36           52.5000   13.4667    Berlin-Ostkreuz                Berlin      Frei
  #00433:      19510101  20260917    48           52.4676   13.4020    Berlin-Tempelhof               Berlin      Frei
  
  # Überschneidung: 1961/01/01 1970/12/31 prüfen
  
  
  #01425:      19480101  19620101   103           50.1214    8.6604    Frankfurt/Main (Feldbergstr.)  Hessen      Frei
  #01421:      19620101  19840901   125           50.1474    8.6746    Frankfurt/Main (Stadt)         Hessen      Frei
  #01420       19810101  20260917   100           50.0259    8.5213    Frankfurt/Main                 Hessen      Frei 
  
  # Überschneidung: 19810101 19840901 prüfen
  
  #01975:      19490101  20260917    11           53.6332    9.9881    Hamburg-Fuhlsbüttel            Hamburg     Frei    
  
  #stations_id von_datum bis_datum Stationshoehe geoBreite geoLaenge   Stationsname                   Bundesland  Abgabe
  #03390       19490101  19920517   529           48.1369   11.7094    München-Riem                   Bayern      Frei
  #03379       19970701  20260917   515           48.1632   11.5429    München-Stadt                  Bayern      Frei
  
  # Keine Überschneidung
  #ABER: 03385 19820101 19990331            515     48.1660   11.5011 München-Nymphenburg                      Bayern                                   Frei   
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

  
  
  
saveRDS(data.d, file = "data.rds")
  
  
  
  

  #####################################################################################################################
  # Data Cleaning for missing values: consider only one NA between two data points in each row and calculate its mean #
  #####################################################################################################################
  
  cat(paste0("Total observed temperature curves (2000-2025): ", dim(data.d)[1],".\n",
               "Before data cleaning: Rows with missing values ",sum(rowSums(is.na(data.d)) > 0), " (",round(sum(rowSums(is.na(data.d)) > 0)/dim(data.d)[1]*100, digits = 2),"%)."))
  
  
  # Berlin:    Total observed temperature curves (2000-2025): 9498. Before data cleaning: Rows with missing values 112 (1.18%).
  # Frankfurt: Total observed temperature curves (2000-2025): 9499. Before data cleaning: Rows with missing values 79  (0.83%).
  # Hamburg:   Total observed temperature curves (2000-2025): 9493. Before data cleaning: Rows with missing values 193 (2.03%).
  # Munich:    Total observed temperature curves (2000-2025): 9499. Before data cleaning: Rows with missing values 92  (0.97%).
  
  
  
  clean.d = data.d |> dplyr::select(4:dim(data.d)[2])
  ind.NA  = data.frame(which(is.na(clean.d), arr.ind = TRUE))
  
  row.neighbor = ind.NA |> group_by(row) |>
                            filter((col - 1) %in% col | (col + 1) %in% col) |>
                              ungroup()   |>
                                select(row) |>
                                  unique()    |>
                                    as.vector()
  
  ind.NA = ind.NA |> filter(!row %in% row.neighbor$row)
  
  # replacing NAs by two neighbors mean
  for(kk in 1:dim(ind.NA)[1]){
      
    if     (ind.NA[kk,2] == 1)  { clean.d[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.d[(ind.NA[kk,1]-1), 144         ]     + clean.d[ ind.NA[kk,1]    , (ind.NA[kk,2]+1)])/2    }
    else if(ind.NA[kk,2] == 144){ clean.d[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.d[ ind.NA[kk,1]   , (ind.NA[kk,2]-1)] + clean.d[(ind.NA[kk,1]+1) , 1           ]    )/2    }
    else                        { clean.d[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.d[ ind.NA[kk,1]   , (ind.NA[kk,2]-1)] + clean.d[ ind.NA[kk,1]    , (ind.NA[kk,2]+1)])/2    }
      
  }
  
  no.neighbor = ind.NA |> group_by(row) |>
                      filter(!(col - 1) %in% col & !(col + 1) %in% col) |>
                        ungroup() |>
                          as.data.frame()
  
  for(kk in 1:dim(no.neighbor)[1]){
    
    if     (ind.NA[kk,2] == 1)  { clean.d[no.neighbor[kk,1], no.neighbor[kk,2]] =  (clean.d[(no.neighbor[kk,1]-1), 144              ]     + clean.d[ no.neighbor[kk,1]    , (no.neighbor[kk,2]+1)])/2    }
    else if(ind.NA[kk,2] == 144){ clean.d[no.neighbor[kk,1], no.neighbor[kk,2]] =  (clean.d[ no.neighbor[kk,1]   , (no.neighbor[kk,2]-1)] + clean.d[(no.neighbor[kk,1]+1) , 1                    ])/2    }
    else                        { clean.d[no.neighbor[kk,1], no.neighbor[kk,2]] =  (clean.d[ no.neighbor[kk,1]   , (no.neighbor[kk,2]-1)] + clean.d[ no.neighbor[kk,1]    , (no.neighbor[kk,2]+1)])/2    }
    
  }
  
  data.d = tibble(data.d[,1:3],clean.d)
  
  
  
  
  
  
  cat(paste0("After data cleaning (one NA, mean of neighbors): Remaining rows with missing values ",sum(rowSums(is.na(data.d)) > 0), " (",round(sum(rowSums(is.na(data.d)) > 0)/dim(data.d)[1]*100, digits = 2),"%)."))  
  
  # Berlin:    After data cleaning (one NA, mean of neighbors): Remaining rows with missing values 88 (0.93%).
  # Frankfurt: After data cleaning (one NA, mean of neighbors): Remaining rows with missing values 51 (0.54%).
  # Hamburg:   After data cleaning (one NA, mean of neighbors): Remaining rows with missing values 98 (1.03%).
  # Munich:    After data cleaning (one NA, mean of neighbors): Remaining rows with missing values 63 (0.66%).
  
  
  
  
  #################################################################################################################################
  # Data Cleaning for missing values: linear interpolation, consider maximum 6 NAs (one hour) between two data points in each row #
  #################################################################################################################################
  
  clean.d = data.d |> dplyr::select(4:dim(data.d)[2])
  ind.NA  = data.frame(which(is.na(clean.d), arr.ind = TRUE))
  
  row.ind = as.numeric(names(table(ind.NA[,1])[table(ind.NA[,1]) <= 6]))  # select each row with maximum 6 NAs
  ind.NA  = ind.NA |> filter(row %in% row.ind) |>  group_by(row)
  
  
  
  # replacing NAs by linear interpolation
  for(kk in unique(ind.NA[,1])$row ){
    
    ind.list = ind.NA |> filter(row == kk)
    ind.min  = min(ind.list$col)
    ind.max  = max(ind.list$col)
    
    if     (ind.min == 1)  { lower = clean.d[(kk-1), 144  ];        upper = clean.d[kk    , (ind.max+1) ] }
    else if(ind.max == 144){ lower = clean.d[ kk   , (k-1)];        upper = clean.d[(kk+1), 1           ] }
    else                   { lower = clean.d[ kk   , (ind.min-1) ]; upper = clean.d[ kk   , (ind.max+1) ] }
    
    x = c(as.numeric(lower), rep(NA, times = (ind.max+1-ind.min) ), as.numeric(upper) )
    
    lin.interpol = round(na_interpolation(x, option = "linear")[-c(1,(ind.max+1-ind.min))], digits = 1)
    clean.d[kk,ind.min:ind.max] =  as.list(lin.interpol) # linear interpolation
    
  }
  
  data.d = tibble(data.d[,1:3],clean.d)
  
  
  
  
  
  cat(paste0("After data cleaning (<= 1h NAs, linear interpolated): Remaining rows with missing values ",sum(rowSums(is.na(data.d)) > 0), " (",round(sum(rowSums(is.na(data.d)) > 0)/dim(data.d)[1]*100, digits = 2),"%)."))  
  
  # Berlin:    After data cleaning (<= 1h NAs, linear interpolated): Remaining rows with missing values 40 (0.42%).
  # Frankfurt: After data cleaning (<= 1h NAs, linear interpolated): Remaining rows with missing values 24 (0.25%).
  # Hamburg:   After data cleaning (<= 1h NAs, linear interpolated): Remaining rows with missing values 50 (0.53%).
  # Munich:    After data cleaning (<= 1h NAs, linear interpolated): Remaining rows with missing values 29 (0.31%).
  
  
  
  
  
  #####################################################################################################################################
  # Data Cleaning for missing values: spline interpolation, consider maximum 18 NAs (three hours) between two data points in each row #
  #####################################################################################################################################
  
  clean.d = data.d |> dplyr::select(4:dim(data.d)[2])
  ind.NA  = data.frame(which(is.na(clean.d), arr.ind = TRUE))
  
  row.ind = as.numeric(names(table(ind.NA[,1])[table(ind.NA[,1]) <= 18]))  # select each row with maximum 6 NAs
  ind.NA  = ind.NA |> filter(row %in% row.ind) |>  group_by(row) |> arrange(row)
  
  # replacing NAs by linear interpolation
  for(kk in unique(ind.NA[,1])$row ){
    
    ind.list = ind.NA |> filter(row == kk)
    ind.min  = min(ind.list$col)
    ind.max  = max(ind.list$col)
    
    if     (ind.min == 1)  { lower = clean.d[(kk-1), 144  ];        upper = clean.d[kk    , (ind.max+1) ] }
    else if(ind.max == 144){ lower = clean.d[ kk   , (k-1)];        upper = clean.d[(kk+1), 1           ] }
    else                   { lower = clean.d[ kk   , (ind.min-1) ]; upper = clean.d[ kk   , (ind.max+1) ] }
    
    x = c(as.numeric(lower), rep(NA, times = (ind.max+1-ind.min) ), as.numeric(upper) )
    
    spline.interpol = round(na_interpolation(x, option = "spline")[-c(1,(ind.max+1-ind.min))], digits = 1)
    clean.d[kk,ind.min:ind.max] =  as.list(spline.interpol) # linear interpolation
    unlist(clean.d[kk,ind.min:ind.max])
  }
  
  data.d = tibble(data.d[,1:3],clean.d)
  
  
  
  
  
  cat(paste0("After data cleaning (1h < and <= 3h NAs, spline interpolated): Rows with missing values ",sum(rowSums(is.na(data.d)) > 0), " (",round(sum(rowSums(is.na(data.d)) > 0)/dim(data.d)[1]*100, digits = 2),"%)."))  
  
  # Berlin:    After data cleaning (1h < and <= 3h NAs, spline interpolated): Rows with missing values 34 (0.36%).
  # Frankfurt: After data cleaning (1h < and <= 3h NAs, spline interpolated): Rows with missing values 13 (0.14%).
  # Hamburg:   After data cleaning (1h < and <= 3h NAs, spline interpolated): Rows with missing values 23 (0.24%).
  # Munich:    After data cleaning (1h < and <= 3h NAs, spline interpolated): Rows with missing values 16 (0.17%).
  
  
  # Removing all over observed curves with NAs
  
  clean.d = data.d |> dplyr::select(4:dim(data.d)[2])
  ind.NA  = data.frame(which(is.na(clean.d), arr.ind = TRUE))
  table(ind.NA[,1])
  
  dense.removed.days = data.d[unique(ind.NA[,1]),1:3]
  
  data.d = data.d[rowSums(is.na(data.d)) == 0,]
  


  ########################################
  # Data Cleaning for false measurements #
  ########################################
  
  clean.d  = data.d |> dplyr::select(4:dim(data.d)[2])
  diff.d   = t(apply(clean.d , 1, diff )) 
  ind      = which(!is.na(diff.d) & abs(diff.d) > 5, arr.ind = TRUE) 
  ind.row  = as.numeric(names(table(ind[,1]))[table(ind[,1]) > 1])
  ind      = ind[ind[,1] %in% ind.row, ]
  ind      = ind[order(ind[, 1]), ]
  
  
  if( k == 1){
    # Berlin: Investigate anomalies manually
  
    # isolated single-point spikes
    unlist(as.vector(clean.d[445  ,c(60:64)]))
    unlist(as.vector(clean.d[1721 ,c(90:94)]))
    unlist(as.vector(clean.d[1728 ,c(103:107)]))
  
    clean.d[445,  62 ] = round((clean.d[445,  61 ]  + clean.d[445,  63 ])/2, digits = 1)
    clean.d[1721, 92 ] = round((clean.d[1721,  91 ] + clean.d[1721,  93])/2, digits = 1)
    clean.d[1728, 105] = round((clean.d[1728, 104 ] + clean.d[1728, 106])/2, digits = 1)
  
    # normal weather
    unlist(as.vector(clean.d[9290 ,c(78:107)]))
  }
  if(k == 2){
    # Frankfurt: Investigate anomalies manually
    
    # isolated single-point spikes
    unlist(as.vector(clean.d[137  ,c(30:33)]))
    unlist(as.vector(clean.d[292  ,c(120:123)]))
    unlist(as.vector(clean.d[2141 ,c(69:73)]))
    unlist(as.vector(clean.d[2350 ,c(20:25)]))
    unlist(as.vector(clean.d[2357 ,c(81:85)]))
    
    
    clean.d[137,  32 ] = round((clean.d[137,  31 ]  + clean.d[137,  33 ])/2, digits = 1)
    clean.d[292, 122 ] = round((clean.d[292, 121 ]  + clean.d[292, 123 ])/2, digits = 1)
    clean.d[2141, 71 ] = round((clean.d[2141, 70 ] + clean.d[2141,  72 ])/2, digits = 1)
    clean.d[2350, 23 ] = round((clean.d[2350, 22 ] + clean.d[2350,  24 ])/2, digits = 1)
    clean.d[2357, 83 ] = round((clean.d[2350, 82 ] + clean.d[2350,  84 ])/2, digits = 1)
    
    # Wrong value sequence
    unlist(as.vector(clean.d[158  ,c(16:25)]))
    unlist(as.vector(clean.d[1960 ,c(1:7)]))
    
    
    clean.d[158,  18:24 ] = as.list(round(na_interpolation(c(unlist(clean.d[158,18]),NA,NA,NA,NA,NA,unlist(clean.d[158,24])) , option = "linear"), digits = 1))
    clean.d[1960, 2 :6 ]  = as.list(round(na_interpolation(c(unlist(clean.d[1960,2]),NA,NA,NA,unlist(clean.d[1960,6])) , option = "linear"), digits = 1))
    
    # normal weather
    unlist(as.vector(clean.d[3848 ,c(108:112)]))
  }
  if(k == 3){
    # Hamburg: Investigate anomalies manually
    
    # isolated single-point spikes
    unlist(as.vector(clean.d[165  ,c(118:122)]))
    unlist(as.vector(clean.d[167  ,c(67:70)]))
    unlist(as.vector(clean.d[193  ,c(30:35)]))
    unlist(as.vector(clean.d[215  ,c(53:58)]))
    unlist(as.vector(clean.d[287  ,c(60:65)]))
    unlist(as.vector(clean.d[490  ,c(4:8)]))
    unlist(as.vector(clean.d[846  ,c(80:85)]))
    unlist(as.vector(clean.d[2696  ,c(40:44)]))
    unlist(as.vector(clean.d[2698  ,c(45:52)]))
    unlist(as.vector(clean.d[2646  ,c(59:63)]))
    unlist(as.vector(clean.d[897  ,c(10:30)]))
    
    clean.d[165, 121 ] = round((clean.d[165, 120 ]  + clean.d[165, 122 ])/2, digits = 1)
    clean.d[167,  69 ] = round((clean.d[167,  68 ]  + clean.d[167,  70 ])/2, digits = 1)
    clean.d[193,  33 ] = round((clean.d[193,  32 ]  + clean.d[193,  34 ])/2, digits = 1)
    clean.d[215,  56 ] = round((clean.d[215,  55 ]  + clean.d[215,  57 ])/2, digits = 1)
    clean.d[287,  63 ] = round((clean.d[287,  62 ]  + clean.d[287,  64 ])/2, digits = 1)
    clean.d[490,  7  ] = round((clean.d[490,   6 ]  + clean.d[490,   8 ])/2, digits = 1)
    clean.d[846,  83 ] = round((clean.d[846,  82 ]  + clean.d[846,  84 ])/2, digits = 1)
    clean.d[2696, 43 ] = round((clean.d[2696, 42 ]  + clean.d[2696, 44 ])/2, digits = 1)
    clean.d[2698, 49 ] = round((clean.d[2698, 48 ]  + clean.d[2698, 50 ])/2, digits = 1)
    clean.d[2646, 61 ] = round((clean.d[2646, 60 ]  + clean.d[2646, 62 ])/2, digits = 1)
    clean.d[897, 16 ]  = round((clean.d[897, 15 ]  + clean.d[897, 17 ])/2, digits = 1)
    clean.d[897, 20 ]  = round((clean.d[897, 19 ]  + clean.d[897, 21 ])/2, digits = 1)
    
    # Wrong value sequence
    unlist(as.vector(clean.d[2698  ,c(50:60)]))
    unlist(as.vector(clean.d[2700  ,c(40:60)]))
    unlist(as.vector(clean.d[2703  ,c(66:71)]))
    unlist(as.vector(clean.d[2695  ,c(25:39)]))
    unlist(as.vector(clean.d[2694  ,c(42:75 )]))
    
    clean.d[2698, 54:56 ] = as.list(c(15.2, 15.8, 15.9))
    clean.d[2700, 47:55 ] = as.list(c(17.0, 17.0, 17.0, 17.0, 17.1, 17.1, 17.4, 17.4, 17.8))
    clean.d[2703, 69 ]    = 27.3
    clean.d[2695, 25:28 ] = as.list(round(na_interpolation(c(unlist(clean.d[2695,25]),NA,NA,unlist(clean.d[2695,28])) , option = "linear"), digits = 1))
    clean.d[2695, 35:38 ] = as.list(round(na_interpolation(c(unlist(clean.d[2695,35]),NA,NA,unlist(clean.d[2695,38])) , option = "linear"), digits = 1))
    clean.d[2694, 42:60 ] = as.list(round(na_interpolation(c(unlist(clean.d[2694,42]),NA,unlist(clean.d[2694,44:47]),NA,unlist(clean.d[2694,49:53]),NA,unlist(clean.d[2694,55:57]),NA,NA,unlist(clean.d[2694,60])) , option = "linear"), digits = 1))
    clean.d[2694, 69:75 ] = as.list(round(na_interpolation(c(unlist(clean.d[2694,69]),NA,unlist(clean.d[2694,71:73]),NA,unlist(clean.d[2694,75])) , option = "linear"), digits = 1))
    
  }
  if(k == 4){
    # Munich: Investigate anomalies manually
    
    # isolated single-point spikes
    unlist(as.vector(clean.d[166,c(138:142)]))
    unlist(as.vector(clean.d[167,c(69:71)]))
    unlist(as.vector(clean.d[174,c(17:19)]))
    unlist(as.vector(clean.d[183,c(53:55)]))
    
    clean.d[166, 141 ] = round((clean.d[166, 140]  + clean.d[166, 142])/2, digits = 1)
    clean.d[167, 70  ] = round((clean.d[166, 69 ]  + clean.d[166, 71 ])/2, digits = 1)
    clean.d[174, 18  ] = round((clean.d[174, 17 ]  + clean.d[174, 19 ])/2, digits = 1)
    clean.d[183, 54  ] = round((clean.d[174, 53 ]  + clean.d[174, 55 ])/2, digits = 1)
    
    
    # Wrong value sequence
    unlist(as.vector(clean.d[787,c(37:50)]))
    unlist(as.vector(clean.d[1532,c(42:60)]))
    
    clean.d[787,  39:45 ] = as.list(round(na_interpolation(c(unlist(clean.d[787,39]),NA,NA,NA,NA,NA,unlist(clean.d[787,45])) , option = "linear"), digits = 1))
    clean.d[1532, 43:52 ] = as.list(round(na_interpolation(c(unlist(clean.d[1532,43]),NA,NA,NA,NA,NA,NA,NA,NA,unlist(clean.d[1532,52])) , option = "linear"), digits = 1))
    
    
    # Remove
    unlist(as.vector(clean.d[1127,c(1:144)]))
    
    clean.d[1127,c(39:66)] = NA 

    # normal weather
  }
  
  
  
  data.d = tibble(data.d[,1:3],clean.d)
  data.d = data.d[rowSums(is.na(data.d)) == 0,] # removed one row in Munich
  
  rm(no.neighbor, x, ind.max, ind.min, kk, lin.interpol, row.ind, ind.row, upper, lower, row.neighbor, ind.NA, ind.list, ind, diff.d, clean.d)
  
  
  
  
  

  
  
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
  
  
  
  
  ##########################################################################################
  # dense data set from 00:00 to 00:00 (145 design points, length 24h)                     #
  ##########################################################################################
  
  data.d.24h = data.d.30h |> select(all_of(keep.cols))
  
  
  ####################################################################################################################
  # dense data set hourly from 00:00 to 00:00 (25 design points, length 24h)                                         #
  # Only usage: Test to assess the cumulative significance of empirical lagged covariance kernel                     #
  # Note: Tensor products of size 145 are computationally too demanding; therefore, hourly aggregated data were used #          
  ####################################################################################################################
  
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
  
  rm(spline.interpol, hour.agg, keep.cols)
  
  
  
  
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
    data.s.part1$MESS_DATUM = data.s.part1$MESS_DATUM - lubridate::hours(1) # MEZ = UTC + 1
    
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
    data.s.part2$MESS_DATUM = data.s.part2$MESS_DATUM - lubridate::hours(1) # MEZ = UTC + 1
    
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
    data.s.part1$MESS_DATUM = data.s.part1$MESS_DATUM - lubridate::hours(1) # MEZ = UTC + 1
    
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
    data.s.part2$MESS_DATUM = data.s.part2$MESS_DATUM - lubridate::hours(1) # MEZ = UTC + 1

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
    data.s$MESS_DATUM = data.s$MESS_DATUM - lubridate::hours(1) # MEZ = UTC + 1
    
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
  
  

  
  print(paste0("Before data cleaning: Remaining rows with missing values ",sum(rowSums(is.na(data.s)) > 0),"."))
  
  # Berlin:    Before data cleaning: Remaining rows with missing values 0.
  # Frankfurt: Before data cleaning: Remaining rows with missing values 1.
  # Hamburg:   Before data cleaning: Remaining rows with missing values 0.
  # Munich:    Before data cleaning: Remaining rows with missing values 0.
  
  
  
  #####################################################################################################################
  # Data Cleaning for missing values: consider only one NA between two data points in each row and calculate its mean #
  #####################################################################################################################
  
  
  
  if(sum(rowSums(is.na(data.s)) > 0) > 0){
    
    clean.s = data.s |> dplyr::select(4:dim(data.s)[2])
    ind.NA  = data.frame(which(is.na(clean.s), arr.ind = TRUE))
  
    row.neighbor = ind.NA |> group_by(row) |>
                    filter((col - 1) %in% col | (col + 1) %in% col) |>
                      ungroup()   |>
                        select(row) |>
                          unique()    |>
                            as.vector()
  
    ind.NA = ind.NA |> filter(!row %in% row.neighbor$row)
  
    # replacing NAs by two neighbors mean
    for(kk in 1:dim(ind.NA)[1]){
    
      if     (ind.NA[kk,2] == 1)  { clean.s[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.s[(ind.NA[kk,1]-1), 24         ]      + clean.s[ ind.NA[kk,1]    , (ind.NA[kk,2]+1)])/2    }
      else if(ind.NA[kk,2] == 24) { clean.s[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.s[ ind.NA[kk,1]   , (ind.NA[kk,2]-1)] + clean.s[(ind.NA[kk,1]+1) , 1           ]    )/2    }
      else                        { clean.s[ind.NA[kk,1], ind.NA[kk,2]] =  (clean.s[ ind.NA[kk,1]   , (ind.NA[kk,2]-1)] + clean.s[ ind.NA[kk,1]    , (ind.NA[kk,2]+1)])/2    }
    
    }
  
    data.s = tibble(data.s[,1:3],clean.s)
  
  
    print(paste0("After data cleaning: Rows with missing values ",sum(rowSums(is.na(data.s)) > 0)))
  
  }
  
  
  # Frankfurt: After data cleaning: Remaining rows with missing values 0.
  
  
  ########################################
  # Data Cleaning for false measurements #
  ########################################
  
  clean.s  = data.s |> dplyr::select(4:dim(data.s)[2])
  diff.s   = t(apply(clean.s , 1, diff )) 
  ind      = which(!is.na(diff.s) & abs(diff.s) > 5, arr.ind = TRUE) 
  ind.row  = as.numeric(names(table(ind[,1]))[table(ind[,1]) > 1])
  ind      = ind[ind[,1] %in% ind.row, ]
  ind      = ind[order(ind[, 1]), ]
  
  
  if( k == 1){
    # Berlin: Investigate anomalies manually
    
    # isolated single-point spikes
    unlist(as.vector(clean.s[7552 ,c(1:6)]))
    clean.s[7552, 2  ] = round((clean.s[7552, 1  ] + clean.s[7552, 3  ])/2, digits = 1)
    
    # normal weather
    unlist(as.vector(clean.s[4289 ,c(7:20)]))
  }
  if(k == 2){
    # Frankfurt: Investigate anomalies manually
    
    # isolated single-point spikes
    unlist(as.vector(clean.s[5461 ,c(2:8)]))
    unlist(as.vector(clean.s[6761 ,c(1:8)]))
    unlist(as.vector(clean.s[6782 ,c(10:18)]))
    
    clean.s[5461, 4 ] = 0.6
    clean.s[6761, 3 ] = round((clean.s[6761, 2 ] + clean.s[6761, 4 ])/2, digits = 1)
    clean.s[6782, 13] = 21.7
    
    #sign error
    unlist(as.vector(clean.s[6536 ,c(10:15)]))
    unlist(as.vector(clean.s[6545 ,c(16:20)]))
    unlist(as.vector(clean.s[6574 ,c(2:5)]))
    
    clean.s[6536, 11:13] = -clean.s[6536, 11:13]
    clean.s[6545, 17 ]   = -clean.s[6545, 17 ]
    clean.s[6574, 3  ]   = -clean.s[6574, 3  ]
    
    # normal weather
    unlist(as.vector(clean.s[3127 ,c(10:19)]))
    unlist(as.vector(clean.s[5339 ,c(12:24)]))
    unlist(as.vector(clean.s[6344 ,c(13:21)]))
    unlist(as.vector(clean.s[6740 ,c(1:15)]))
    
  }
  if(k == 3){
    # Hamburg: Investigate anomalies manually

    # Wrong value sequence
    unlist(as.vector(clean.s[3157 ,c(10:15)]))
    unlist(as.vector(clean.s[4169 ,c(1:10)]))
    
    clean.s[3157, 12 ]   = 18.0
    clean.s[4169,  2 ]   = 9.4
    
    #sign error
    unlist(as.vector(clean.s[1121 ,c(10:15)]))
    clean.s[1121, 13 ]   = -clean.s[1121, 13 ]
  
    # normal weather
    unlist(as.vector(clean.s[2010 ,c(16:20)]))
    unlist(as.vector(clean.s[2804 ,c(5:10)]))
    unlist(as.vector(clean.s[3821 ,c(13:20)]))
    unlist(as.vector(clean.s[4230 ,c(10:17)]))
    unlist(as.vector(clean.s[6037 ,c(6:10)]))
    unlist(as.vector(clean.s[7463 ,c(15:20)]))
    
  }
  if(k == 4){
    # Munich: Investigate anomalies manually
    
    # normal weather
    unlist(as.vector(clean.s[6033 ,c(5:17)]))
    unlist(as.vector(clean.s[5167 ,c(9:14)]))
    unlist(as.vector(clean.s[5037 ,c(9:13)]))
    unlist(as.vector(clean.s[4347 ,c(8:17)]))
    unlist(as.vector(clean.s[2578 ,c(19:23)]))
    unlist(as.vector(clean.s[1241 ,c(5:11)]))
    unlist(as.vector(clean.s[195 ,c(11:17)]))
    unlist(as.vector(clean.s[4583 ,c(11:20)]))

  }
  
  
  
  
  
  rm(x, ind.max, ind.min, kk, lin.interpol, row.ind, ind.row, upper, lower, row.neighbor, ind.NA, ind.list, ind, diff.s, clean.s)
  
  
  
  
  
  

  
  
  
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


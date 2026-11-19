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


# constant difference

for(j in 1:12){
  
  i = 3
  
  delta_constant = TRUE
  scenario3      = FALSE    
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
  i = 2
  
  delta_constant = TRUE
  scenario3      = FALSE    
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
  
  i = 1
  
  delta_constant = TRUE
  scenario3      = FALSE     
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
}



# non-constant difference

for(j in 1:12){

    i = 3

    delta_constant = FALSE
    scenario3      = FALSE    
    integral       = TRUE   
  
    source("Simulation/Comparison/both_estimations.R")
    
    i = 2

    delta_constant = FALSE
    scenario3      = FALSE    
    integral       = TRUE   
  
    source("Simulation/Comparison/both_estimations.R")

    
    i = 1

    delta_constant = FALSE
    scenario3      = FALSE     
    integral       = TRUE   
  
    source("Simulation/Comparison/both_estimations.R")
    
}





# non-constant difference with more features than mean functions

for(j in 1:12){
  
  i = 3
  
  delta_constant = FALSE
  scenario3      = TRUE    
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
  i = 2
  
  delta_constant = FALSE
  scenario3      = TRUE    
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
  
  i = 1
  
  delta_constant = FALSE
  scenario3      = TRUE     
  integral       = TRUE   
  
  source("Simulation/Comparison/both_estimations.R")
  
}

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


# not centered difference function

for(j in 2:11){
  
  i = 1
  
  
  delta_constant = TRUE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
}

for(j in 2:11){
  
  i = 2
  
  
  delta_constant = TRUE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
}



for(j in 2:11){
  
  i = 3
  
  
  delta_constant = TRUE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = FALSE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  

}






# centered difference function

for(j in 2:11){
  
  i = 1
  
  
  delta_constant = TRUE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
}

for(j in 2:11){
  
  i = 2
  
  
  delta_constant = TRUE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
}



for(j in 2:11){
  
  i = 3
  
  
  delta_constant = TRUE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
  
  delta_constant = FALSE    
  integral       = TRUE   
  source("Simulation/Coverage and power/Coverage_and_power_of_Test.R")
  
}
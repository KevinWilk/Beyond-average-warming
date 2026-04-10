#!/bin/bash
#SBATCH --job-name=Simulation_future
#SBATCH --partition=mqtest
#SBATCH --time=9-23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=245000M
#SBATCH --output=Logs/Simulation_future_%j.out

module load deps/other/gcc/14.2.0
module load R/4.1.2

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
Rscript "weather temperature/confidence bands/constucting_confidence_bands.R" 
#Rscript "weather temperature/bandwidth selection/bandwidth_selection.R" 
#Rscript "weather temperature/kernel weights/calculation_weights.R" 
#Rscript "weather temperature/long run kernel/lagged_covariance.R" 

#Rscript "Simulation/Coverage and power/Run_simulation.R" 
#Rscript "Simulation/Coverage and power/Coverage_and_power_of_Test.R"
#Rscript "Simulation/optimal bandwidth/LagCov_bandwidth_selection.R" 
#Rscript "Simulation/optimal bandwidth/Mean_bandwidth_selection.R" 
#Rscript "Simulation/weights/Calculating _weights.R" 

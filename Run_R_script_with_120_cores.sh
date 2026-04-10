#!/bin/bash
#SBATCH --job-name=Simulation_future
#SBATCH --partition=mqtest
#SBATCH --time=2-23:59:00
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
Rscript Run_simulation.R 

#!/bin/bash
#SBATCH -J int
#SBATCH --partition=patralab
#SBATCH --time=7-0:00:00 
#SBATCH -n 10
#SBATCH --mem=500Gb
#SBATCH --output=batch_rscript_%j.out #saving standard output to file
#SBATCH --error=batch_rscript_%j.err #saving standard error to file
 
module purge
module load R/4.0.0

Rscript --no-save $1

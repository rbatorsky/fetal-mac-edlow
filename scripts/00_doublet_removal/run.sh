#!/bin/bash
#SBATCH -J doublet_removal
#SBATCH --partition=patralab
#SBATCH --time=1-0:00:00 
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64Gb
#SBATCH --output=doublet_remove_obsctr_%j.out #saving standard output to file
#SBATCH --error=doublet_remove_obsctr_%j.err #saving standard error to file
 
module purge
module load gcc/7.3.0
module load anaconda/3
module load R/4.0.0

sample_name=$1
cellranger_directory=/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/all_cellranger/
output_path=/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/analysis/doublet_removal/

Rscript --no-save 00_doublet_removal.R $sample_name $cellranger_directory $output_path


#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem 20G
#SBATCH -t 0-5:00
#SBATCH -o slurm/MCCdna-%j.out
#SBATCH -e slurm/MCCdna-%j.err

module load matlab

# change filename to Rmd to be knitted. 
# Make sure ncores in Rmd matches -n batch param above
R -e "rmarkdown::render('20191121-analysis.Rmd')"

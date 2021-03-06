#!/bin/bash
#SBATCH -J MCCrna
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem 20G
#SBATCH -t 0-4:00
#SBATCH -o _slurm/MCCrna-%j.out
#SBATCH -e _slurm/MCCrna-%j.err

#export RSTUDIO_PANDOC="/n/helmod/apps/centos7/Core/rstudio/1.1.453-fasrc01/bin/pandoc"

# change filename to Rmd to be knitted. 
# Make sure ncores in Rmd matches -n batch param above
R -e "rmarkdown::render('20191111-analysis.Rmd')"

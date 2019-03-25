#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 20G
#SBATCH -t 0-1:00
#SBATCH -o slurm/MCCdna-%j.out
#SBATCH -e slurm/MCCdna-%j.err

export RSTUDIO_PANDOC="/n/helmod/apps/centos7/Core/rstudio/1.1.453-fasrc01/bin/pandoc"

# change filename to Rmd to be knitted. 
# Make sure ncores in Rmd matches -n batch param above
R -e "rmarkdown::render('20190325-analysis.Rmd')"

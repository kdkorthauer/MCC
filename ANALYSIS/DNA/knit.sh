#!/bin/bash
#PBS -l walltime=06:00:00,select=1:ncpus=10:mem=20gb
#PBS -N MCCdna
#PBS -A st-kdkortha-1
#PBS -o /scratch/st-kdkortha-1/_pbs/MCC/MCC_DNA.out
#PBS -e /scratch/st-kdkortha-1/_pbs/MCC/MCC_DNA.err


export RSTUDIO_PANDOC="/home/kdkortha/bin/rstudio-1.2.5033/bin/pandoc"
module load matlab

# change filename to Rmd to be knitted. 
# Make sure ncores in Rmd matches -n batch param above
R -e "rmarkdown::render('20191121-analysis.Rmd')"

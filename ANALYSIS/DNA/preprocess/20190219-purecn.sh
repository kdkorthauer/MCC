#!/bin/bash
#SBATCH -J purecn
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 70G
#SBATCH -t 0-18:00

R CMD BATCH --quiet --no-restore --no-save 20190219-purecn.R ../slurm/20190219-purecn.Rout
#!/bin/bash
#SBATCH -J purecn
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 70G
#SBATCH -t 0-20:00

#SBATCH -o ../slurm/purecn_%j.out    # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../slurm/purecn_%j.err     # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-user=kdkorthauer@gmail.com

R CMD BATCH --quiet --no-restore --no-save 20190219-purecn.R ../slurm/20190219-purecn.Rout